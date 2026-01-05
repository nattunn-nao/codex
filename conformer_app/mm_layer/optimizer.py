from __future__ import annotations

"""
mm_layer.optimizer

- builder.build_clustered_multiconformer_mol() を用いて
  「フラグメントを shell 上に配置したクラスター構造」を複数生成
- MMFF / UFF 等で MM 最適化
- エネルギー重複除去 & エネルギー窓フィルタ & 上位 max_keep 選抜
- XYZ 出力と簡易ログ出力を行う

公開関数:
    run_mm_for_job(mol_id, smiles, embed_cfg, mm_cfg, out_cfg) -> MMResult
"""

from dataclasses import dataclass
from typing import List, Tuple

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from conformer_app.core.config import EmbedConfig, MMConfig, OutputConfig
from conformer_app.app_io.paths import id_root_dir, mm_dir
from conformer_app.app_io.xyz_io import write_single_xyz
from .builder import build_clustered_multiconformer_mol


# ------------------------------
# 結果コンテナ
# ------------------------------

@dataclass
class MMResult:
    mol_id: str
    smiles: str
    n_generated: int = 0          # 生成したクラスター構造数
    n_optimized: int = 0          # MM 最適化に成功した個数
    n_after_dedup: int = 0        # エネルギー重複除去 & 窓後の個数
    n_kept: int = 0               # 実際に XYZ に書き出した個数
    mm_xyz_files: List[Path] | None = None  # UMA に渡す XYZ ファイル一覧


# ------------------------------
# 内部ユーティリティ
# ------------------------------

def _build_forcefield(mol: Chem.Mol, conf_id: int, ff_mode: str, max_iters: int):
    """
    RDKit 力場オブジェクトを構築する。
    ff_mode:
        - "MMFF94s->UFF"
        - "MMFF94"
        - "UFF"
    """
    ff_mode = (ff_mode or "MMFF94s->UFF").upper()

    if ff_mode.startswith("MMFF94"):
        mmff_props = AllChem.MMFFGetMoleculeProperties(
            mol,
            mmffVariant="MMFF94S" if "S" in ff_mode else "MMFF94",
        )
        if mmff_props is not None:
            ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)
            if ff is not None:
                ff.Initialize()
                ff.SetMaxIterations(int(max_iters))
                return ff
        # フォールバック UFF
        if "UFF" in ff_mode:
            uff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
            if uff is not None:
                uff.Initialize()
                uff.SetMaxIterations(int(max_iters))
                return uff
        return None

    # UFF 指定
    if ff_mode == "UFF":
        uff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        if uff is not None:
            uff.Initialize()
            uff.SetMaxIterations(int(max_iters))
            return uff

    # デフォルト: MMFF94s->UFF
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94S")
    if mmff_props is not None:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)
        if ff is not None:
            ff.Initialize()
            ff.SetMaxIterations(int(max_iters))
            return ff

    uff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
    if uff is not None:
        uff.Initialize()
        uff.SetMaxIterations(int(max_iters))
    return uff


def _mm_optimize_all_confs(
    mol: Chem.Mol,
    conf_ids: List[int],
    mm_cfg: MMConfig,
) -> List[Tuple[int, float]]:
    """
    全ての Conformer について MM 最適化を行い、
    (conf_id, energy_kcal) のリストを返す（失敗したものは除外）。
    """
    results: List[Tuple[int, float]] = []
    for cid in conf_ids:
        ff = _build_forcefield(mol, cid, mm_cfg.ff_mode, mm_cfg.max_iters)
        if ff is None:
            continue
        try:
            _ = ff.Minimize()
            e = float(ff.CalcEnergy())  # RDKit 力場のエネルギー（kcal/mol 相当）
            results.append((cid, e))
        except Exception:
            continue
    return results


def _energy_dedup_and_window(
    energies: List[Tuple[int, float]],
    equal_tol: float,
    energy_window: float,
) -> List[Tuple[int, float]]:
    """
    - エネルギー昇順にソート
    - ΔE ≤ equal_tol を同値とみなして代表 1 つだけ残す
    - Emin からのエネルギー窓 (E ≤ Emin + window) を適用
    """
    if not energies:
        return []

    energies_sorted = sorted(energies, key=lambda x: x[1])
    kept: List[Tuple[int, float]] = []

    base_conf_id, base_e = energies_sorted[0]
    kept.append((base_conf_id, base_e))

    for cid, e in energies_sorted[1:]:
        if e - base_e > energy_window:
            # 窓から外れたら打ち切り
            break
        last_e = kept[-1][1]
        if abs(e - last_e) <= equal_tol:
            continue
        kept.append((cid, e))

    return kept


def _center_of_mass(mol: Chem.Mol, conf_id: int) -> tuple[float, float, float]:
    """コンフォマーごとの重心座標を返す。"""
    conf = mol.GetConformer(conf_id)
    n = mol.GetNumAtoms()
    xs = ys = zs = 0.0
    for i in range(n):
        pos = conf.GetAtomPosition(i)
        xs += float(pos.x)
        ys += float(pos.y)
        zs += float(pos.z)
    return xs / n, ys / n, zs / n


def _translate_to_com_origin(mol: Chem.Mol, conf_id: int) -> None:
    """重心が原点に来るように座標を平行移動する。"""
    cx, cy, cz = _center_of_mass(mol, conf_id)
    conf = mol.GetConformer(conf_id)
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        pos.x -= cx
        pos.y -= cy
        pos.z -= cz
        conf.SetAtomPosition(i, pos)


def _write_selected_xyz(
    mol: Chem.Mol,
    mol_id: str,
    selected: List[Tuple[int, float]],
    out_dir: Path,
    move_to_com_origin: bool,
) -> List[Path]:
    """
    選抜された Conformer 群を XYZ として out_dir に書き出す。
    ファイル名: confNNN.xyz（energy をコメントに付与）
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    for idx, (cid, e) in enumerate(selected, start=1):
        if move_to_com_origin:
            _translate_to_com_origin(mol, cid)

        conf = mol.GetConformer(cid)
        coords = []
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append(
                (atom.GetSymbol(), float(pos.x), float(pos.y), float(pos.z))
            )

        fname = out_dir / f"conf{idx:03d}.xyz"
        comment = f"{mol_id} conf={idx} energy={e:.6f} kcal/mol"
        write_single_xyz(fname, coords, comment=comment)
        written.append(fname)

    return written


def _append_global_mm_log(
    out_cfg: OutputConfig,
    mol_id: str,
    n_generated: int,
    n_optimized: int,
    n_after_dedup: int,
    n_kept: int,
) -> None:
    """
    ./<out_base>/logs/mm_summary.log に 1 行追記する。
    """
    base = Path(out_cfg.base_dir).resolve()
    logs_dir = base / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / "mm_summary.log"

    line = (
        f"{mol_id}\t"
        f"generated={n_generated}\t"
        f"optimized={n_optimized}\t"
        f"after_dedup={n_after_dedup}\t"
        f"kept={n_kept}\n"
    )
    try:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(line)
    except Exception:
        # ログ失敗は致命的ではないので握りつぶす
        pass


# ------------------------------
# 公開関数
# ------------------------------

def run_mm_for_job(
    mol_id: str,
    smiles: str,
    embed_cfg: EmbedConfig,
    mm_cfg: MMConfig,
    out_cfg: OutputConfig,
) -> MMResult:
    """
    1 つの ID (= mol_id) について:
      1. builder.build_clustered_multiconformer_mol() でクラスター構造群を生成
      2. MM 力場で最適化
      3. エネルギー昇順ソート＋重複除去＋エネルギー窓
      4. 上位 max_keep を XYZ に出力
    """
    result = MMResult(mol_id=mol_id, smiles=smiles, mm_xyz_files=[])

    # ルートディレクトリ
    id_dir = id_root_dir(mol_id, out_cfg)
    mm_out_dir = mm_dir(mol_id, out_cfg)
    id_dir.mkdir(parents=True, exist_ok=True)
    mm_out_dir.mkdir(parents=True, exist_ok=True)

    # 1) クラスター構造生成（フラグメントシェル配置＋原子間距離チェック）
    mol, conf_ids = build_clustered_multiconformer_mol(smiles, embed_cfg)
    if not conf_ids:
        raise RuntimeError(f"[{mol_id}] クラスター構造の生成に失敗しました。")
    result.n_generated = len(conf_ids)

    # 2) MM 最適化
    energies = _mm_optimize_all_confs(mol, conf_ids, mm_cfg)
    if not energies:
        raise RuntimeError(f"[{mol_id}] MM 最適化に全て失敗しました。")
    result.n_optimized = len(energies)

    # 3) エネルギー重複除去 & 窓
    equal_tol = float(getattr(mm_cfg, "equal_tol", 0.003))
    energy_window = float(getattr(mm_cfg, "energy_window", 10.0))
    deduped = _energy_dedup_and_window(
        energies,
        equal_tol=equal_tol,
        energy_window=energy_window,
    )
    if not deduped:
        raise RuntimeError(f"[{mol_id}] エネルギー重複除去・窓適用後に構造が残りませんでした。")
    result.n_after_dedup = len(deduped)

    # 4) 上位 max_keep を選抜して XYZ 出力
    max_keep = int(getattr(mm_cfg, "max_keep", 10))
    selected = deduped[:max_keep]
    move_to_com_origin = bool(out_cfg.move_to_com_origin)

    written = _write_selected_xyz(
        mol,
        mol_id=mol_id,
        selected=selected,
        out_dir=mm_out_dir,
        move_to_com_origin=move_to_com_origin,
    )
    result.n_kept = len(written)
    result.mm_xyz_files = written

    # グローバル MM ログに追記
    _append_global_mm_log(
        out_cfg,
        mol_id=mol_id,
        n_generated=result.n_generated,
        n_optimized=result.n_optimized,
        n_after_dedup=result.n_after_dedup,
        n_kept=result.n_kept,
    )

    return result
