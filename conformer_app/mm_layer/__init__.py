from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import math
import random
import shutil

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Geometry import Point3D

from conformer_app.core.config import EmbedConfig, MMConfig, OutputConfig
from conformer_app.app_io.paths import (
    id_root_dir,
    mm_dir,
    post_mm_dir,
    logs_dir,
)

# フラグメント同士の「表面間距離」の目安 [Å]
_FRAG_GAP_MIN = 0.5   # 近すぎない下限
_FRAG_GAP_MAX = 1.5   # 離れすぎない上限


@dataclass
class MMResult:
    """1 ID の MM ステージの簡易サマリ"""
    n_generated: int          # 生成した「クラスター構造」（= initial conformers）の数
    n_optimized: int          # MM 最適化に成功したコンフォマーの数
    xyz_files: List[str]      # UMA に渡す XYZ（post_MM）の一覧


# =========================
# 基本ユーティリティ
# =========================

def _split_smiles(smiles: str) -> List[str]:
    """'.' で SMILES をフラグメントに分割"""
    parts = [p for p in smiles.replace(" ", "").split(".") if p]
    if not parts:
        raise ValueError(f"フラグメントが見つかりません: {smiles}")
    return parts


def _sanitize_smiles_to_mol(smiles: str) -> Chem.Mol:
    """
    単一フラグメント SMILES → H 付き・sanitize 済み Mol
    （多フラグメントは別ルートで処理する）
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES から Mol を生成できませんでした: {smiles}")
    mol = Chem.AddHs(mol)
    if not mol.GetNumAtoms():
        raise ValueError(f"H 付加後の分子が空です: {smiles}")

    # 明示的にサニタイズして RingInfo や各種プロパティを初期化
    Chem.SanitizeMol(mol)
    mol.UpdatePropertyCache(strict=False)

    return mol


def _make_etkdg_params(embed_cfg: EmbedConfig) -> AllChem.EmbedParameters:
    """EmbedConfig から ETKDG のパラメータを構築"""
    ver = (embed_cfg.etkdg_version or "ETKDGv3").upper()
    if ver == "ETKDGV3":
        params = AllChem.ETKDGv3()
    elif ver == "ETKDGV2":
        params = AllChem.ETKDGv2()
    else:
        params = AllChem.ETKDG()

    params.pruneRmsThresh = (
        float(embed_cfg.prune_rms) if embed_cfg.prune_rms > 0 else -1.0
    )
    params.randomSeed = int(embed_cfg.random_seed)
    params.useSmallRingTorsions = bool(embed_cfg.small_ring)
    params.useMacrocycleTorsions = bool(embed_cfg.macrocycle)
    return params


# =========================
# フラグメント用: コンフォマー生成 & シェル配置
# =========================

def _compute_radius(coords: np.ndarray) -> float:
    """
    coords: (N, 3), すでに重心を原点に移した座標を想定。
    原点から最も遠い点までの距離を「半径」とみなす。
    """
    if coords.size == 0:
        return 0.0
    r2 = np.sum(coords * coords, axis=1)
    return float(math.sqrt(float(r2.max())))


def _random_unit_vector() -> Tuple[float, float, float]:
    """3次元の一様乱数方向ベクトル"""
    z = random.uniform(-1.0, 1.0)
    t = random.uniform(0.0, 2.0 * math.pi)
    r_xy = math.sqrt(max(0.0, 1.0 - z * z))
    x = r_xy * math.cos(t)
    y = r_xy * math.sin(t)
    return (x, y, z)


def _dist(p: Tuple[float, float, float], q: Tuple[float, float, float]) -> float:
    """2点間距離"""
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    dz = p[2] - q[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _build_fragment_conformers(
    smiles: str,
    n_confs: int,
    params: AllChem.EmbedParameters,
) -> Tuple[Chem.Mol, List[Dict[str, np.ndarray]], int]:
    """
    1つのフラグメント SMILES に対して、
    最大 n_confs 個のコンフォマーを生成し、
    重心を原点に平行移動した座標と半径を返す。

    戻り値:
        frag_mol : Chem.Mol (H 付き)
        conf_list : list of {"coords": (N,3) np.ndarray, "radius": float}
        heavy : 重原子数
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES 読み込み失敗: {smiles}")
    mol = Chem.AddHs(mol)

    # H 追加後のフラグメントもここでサニタイズして RingInfo 等を初期化
    Chem.SanitizeMol(mol)
    mol.UpdatePropertyCache(strict=False)

    heavy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)

    conf_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=int(n_confs),
        params=params,
    )
    conf_ids = list(conf_ids)
    if not conf_ids:
        raise RuntimeError(f"コンフォマー埋め込みに失敗: {smiles}")

    # ざっくり UFF で形だけ整える
    AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0)

    conf_list: List[Dict[str, np.ndarray]] = []

    for cid in conf_ids:
        conf = mol.GetConformer(cid)
        n_atoms = mol.GetNumAtoms()
        coords = np.zeros((n_atoms, 3), dtype=float)
        for i in range(n_atoms):
            pos = conf.GetAtomPosition(i)
            coords[i] = [pos.x, pos.y, pos.z]

        # 重心を原点に移動（質量は無視して単純平均）
        com = coords.mean(axis=0)
        centered = coords - com
        radius = _compute_radius(centered)
        conf_list.append({"coords": centered, "radius": radius})

    return mol, conf_list, heavy


def _place_fragments_on_shell(
    radii: List[float],
    gap_min: float,
    gap_max: float,
) -> List[Tuple[float, float, float]]:
    """
    radii: [r_main, r_1, r_2, ...] （0番目が代表フラグメント）
    代表フラグメントを原点に置き、残りを shell 内に配置する。
    """
    if gap_max <= gap_min:
        raise ValueError("gap_max は gap_min より大きくしてください。")

    n = len(radii)
    if n == 0:
        return []

    r_main = float(radii[0])
    R_inner = r_main + gap_min
    R_outer = r_main + gap_max

    placements: List[Tuple[float, float, float]] = [(0.0, 0.0, 0.0)]  # 代表は原点

    if n == 1:
        return placements

    max_trials = 300

    for i in range(1, n):
        r_i = float(radii[i])
        placed = False

        for _ in range(max_trials):
            # shell 内の「最近接末端半径」をランダムに選ぶ
            r_surface = random.uniform(R_inner, R_outer)
            # 中心距離 = 末端半径 + 自身の半径
            d = r_surface + r_i

            ux, uy, uz = _random_unit_vector()
            cand = (d * ux, d * uy, d * uz)

            ok = True
            # 他フラグメントとの重なりをチェック
            for j in range(i):
                r_j = float(radii[j])
                d_ij = _dist(cand, placements[j])
                d_min = r_i + r_j + gap_min
                if d_ij < d_min:
                    ok = False
                    break

            if ok:
                placements.append(cand)
                placed = True
                break

        if not placed:
            # 条件が厳しすぎたら外側に逃がす
            ux, uy, uz = _random_unit_vector()
            d = R_outer + r_i + gap_min
            placements.append((d * ux, d * uy, d * uz))

    return placements


# =========================
# 初期コンフォマー生成（単一 / 多フラグメント両対応）
# =========================

def _prepare_initial_confs_for_smiles(
    smiles: str,
    embed_cfg: EmbedConfig,
) -> Tuple[Chem.Mol, List[int]]:
    """
    単一フラグメント:
        → 従来どおり ETKDG で multi-Conf

    多フラグメント ('.' を含む):
        → フラグメントごとにコンフォマーをプールし、
           ランダムに組み合わせて「クラスター構造」を embed_cfg.num_confs 個生成。
    """
    frag_smis = _split_smiles(smiles)

    # 再現性が欲しい場合用（UI から random_seed が渡ってくる）
    if embed_cfg.random_seed >= 0:
        random.seed(int(embed_cfg.random_seed))

    # --- 単一フラグメント: 従来型 ---
    if len(frag_smis) == 1:
        mol = _sanitize_smiles_to_mol(smiles)
        params = _make_etkdg_params(embed_cfg)
        conf_ids = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=int(embed_cfg.num_confs),
            params=params,
        )
        conf_ids = list(conf_ids)
        if not conf_ids:
            raise RuntimeError("ETKDG によるコンフォマー生成に失敗しました。")
        return mol, conf_ids

    # --- 多フラグメント: クラスター生成アルゴリズム ---
    params = _make_etkdg_params(embed_cfg)

    frag_infos = []
    frag_mols: List[Chem.Mol] = []

    for fs in frag_smis:
        frag_mol, conf_list, heavy = _build_fragment_conformers(
            fs,
            n_confs=int(embed_cfg.num_confs),
            params=params,
        )
        frag_infos.append(
            {
                "smiles": fs,
                "mol": frag_mol,
                "confs": conf_list,   # list of {"coords", "radius"}
                "heavy": heavy,
            }
        )
        frag_mols.append(frag_mol)

    # CombineMols で「結合の無いフラグメント集合」として 1 Mol に統合
    combined = frag_mols[0]
    for fm in frag_mols[1:]:
        combined = Chem.CombineMols(combined, fm)

    # CombineMols の戻り値はサニタイズされておらず RingInfo も未初期化なので、
    # ここでサニタイズして RingInfo / プロパティキャッシュを初期化しておく。
    rw_combined = Chem.RWMol(combined)
    Chem.SanitizeMol(rw_combined)
    combined = rw_combined.GetMol()
    combined.UpdatePropertyCache(strict=False)

    total_atoms = combined.GetNumAtoms()
    conf_ids: List[int] = []

    n_structures = int(embed_cfg.num_confs)

    for _ in range(n_structures):
        # 各フラグメントからコンフォマーをランダムに1つ選択
        chosen = []
        for info in frag_infos:
            conf_data = random.choice(info["confs"])
            chosen.append(
                {
                    "coords": conf_data["coords"],   # (n_i,3)
                    "radius": conf_data["radius"],
                }
            )

        # 一番重いフラグメントを代表構造に
        main_idx = max(range(len(frag_infos)), key=lambda i: frag_infos[i]["heavy"])
        ordered_indices = [main_idx] + [
            i for i in range(len(frag_infos)) if i != main_idx
        ]

        radii = [float(chosen[i]["radius"]) for i in ordered_indices]
        offsets = _place_fragments_on_shell(
            radii,
            gap_min=_FRAG_GAP_MIN,
            gap_max=_FRAG_GAP_MAX,
        )

        # CombineMols は「フラグメントをそのまま連結」するので、
        # atom の並び順は [frag0 の全原子, frag1 の全原子, ...] の順になる。
        coords_all = np.zeros((total_atoms, 3), dtype=float)
        atom_cursor = 0
        for order_pos, frag_idx in enumerate(ordered_indices):
            conf_data = chosen[frag_idx]
            base = np.array(conf_data["coords"], dtype=float)
            off = np.array(offsets[order_pos], dtype=float).reshape(1, 3)
            shifted = base + off
            n_i = shifted.shape[0]
            coords_all[atom_cursor:atom_cursor + n_i, :] = shifted
            atom_cursor += n_i

        # RDKit Conformer として登録
        conf = Chem.Conformer(total_atoms)
        for i in range(total_atoms):
            x, y, z = coords_all[i]
            conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
        cid = combined.AddConformer(conf, assignId=True)
        conf_ids.append(cid)

    return combined, conf_ids


# =========================
# MM 最適化 & 重複除去
# =========================

def _optimize_one_conf_mmff(
    mol: Chem.Mol,
    conf_id: int,
    max_iters: int,
    variant: str = "MMFF94s",
) -> float:
    """MMFF で 1 コンフォマーを最適化し、エネルギーを返す。"""
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=variant)
    if props is None:
        raise RuntimeError("MMFF のパラメータを取得できませんでした。")

    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    if ff is None:
        raise RuntimeError("MMFF のフォースフィールド生成に失敗しました。")

    ff.Initialize()
    ff.Minimize(maxIts=int(max_iters))
    e = float(ff.CalcEnergy())
    return e


def _optimize_one_conf_uff(
    mol: Chem.Mol,
    conf_id: int,
    max_iters: int,
) -> float:
    """UFF で 1 コンフォマーを最適化し、エネルギーを返す。"""
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
    if ff is None:
        raise RuntimeError("UFF のフォースフィールド生成に失敗しました。")

    ff.Initialize()
    ff.Minimize(maxIts=int(max_iters))
    e = float(ff.CalcEnergy())
    return e


def _mm_optimize_all_confs(
    mol: Chem.Mol,
    conf_ids: List[int],
    cfg: MMConfig,
) -> Dict[int, float]:
    """全コンフォマーを MM 最適化し、conf_id → エネルギーの辞書を返す。"""
    energies: Dict[int, float] = {}

    for cid in conf_ids:
        if cfg.ff_mode == "MMFF94s->UFF":
            try:
                e = _optimize_one_conf_mmff(mol, cid, cfg.max_iters, variant="MMFF94s")
            except Exception:
                e = _optimize_one_conf_uff(mol, cid, cfg.max_iters)
        elif cfg.ff_mode == "MMFF94":
            e = _optimize_one_conf_mmff(mol, cid, cfg.max_iters, variant="MMFF94")
        elif cfg.ff_mode == "UFF":
            e = _optimize_one_conf_uff(mol, cid, cfg.max_iters)
        else:
            raise ValueError(f"未知の ff_mode: {cfg.ff_mode}")

        energies[cid] = float(e)

    return energies


def _calc_rmsd_between_confs(mol: Chem.Mol, cid_ref: int, cid_prb: int) -> float:
    """2 コンフォマー間 RMSD（最適重ね合わせ）"""
    try:
        rmsd = rdMolAlign.GetBestRMS(mol, mol, refId=cid_ref, prbId=cid_prb)
        return float(rmsd)
    except Exception:
        # フォールバック: 単純な座標 RMSD（重ね合わせなし）
        conf_ref = mol.GetConformer(cid_ref)
        conf_prb = mol.GetConformer(cid_prb)
        n = mol.GetNumAtoms()
        coords_ref = np.zeros((n, 3), dtype=float)
        coords_prb = np.zeros((n, 3), dtype=float)
        for i in range(n):
            p1 = conf_ref.GetAtomPosition(i)
            p2 = conf_prb.GetAtomPosition(i)
            coords_ref[i] = [p1.x, p1.y, p1.z]
            coords_prb[i] = [p2.x, p2.y, p2.z]
        diff = coords_ref - coords_prb
        return float(np.sqrt((diff * diff).sum() / n))


def _select_unique_confs(
    mol: Chem.Mol,
    energies: Dict[int, float],
    sorted_items: List[Tuple[int, float]],
    cfg: MMConfig,
) -> List[int]:
    """
    エネルギー & RMSD 重複除去を行い、選抜するコンフォマーの conf_id リストを返す。

    - equal_tol (kcal/mol) 以内で
    - RMSD < rmsd_min のものは重複とみなして捨てる
    - 全体として E0 + energy_window を超えるものは見ない
    - max_keep までを選抜
    """
    if not sorted_items:
        return []

    equal_tol = float(getattr(cfg, "equal_tol", 0.003))
    energy_window = float(getattr(cfg, "energy_window", 10.0))
    rmsd_min = float(getattr(cfg, "rmsd_min", 0.8))
    max_keep = int(getattr(cfg, "max_keep", 10))

    kept: List[int] = []

    cid0, E0 = sorted_items[0]
    kept.append(cid0)

    for cid, E in sorted_items[1:]:
        if E - E0 > energy_window:
            break

        is_dup = False
        for kcid in kept:
            dE = abs(E - energies[kcid])
            if dE > equal_tol:
                continue
            rmsd = _calc_rmsd_between_confs(mol, kcid, cid)
            if rmsd < rmsd_min:
                is_dup = True
                break

        if is_dup:
            continue

        kept.append(cid)
        if len(kept) >= max_keep:
            break

    return kept


def _move_positions_to_com(mol: Chem.Mol, conf_id: int) -> np.ndarray:
    """
    指定コンフォマーの座標を取得し、重心が原点にくるように平行移動した座標配列を返す。
    （mol 自体の座標は書き換えない）
    """
    conf = mol.GetConformer(conf_id)
    n = mol.GetNumAtoms()
    coords = np.zeros((n, 3), dtype=float)
    masses = np.zeros((n,), dtype=float)

    pt = Chem.GetPeriodicTable()

    for i in range(n):
        pos = conf.GetAtomPosition(i)
        coords[i, 0] = pos.x
        coords[i, 1] = pos.y
        coords[i, 2] = pos.z
        masses[i] = pt.GetAtomicWeight(mol.GetAtomWithIdx(i).GetAtomicNum())

    total_mass = masses.sum()
    if total_mass <= 0:
        com = coords.mean(axis=0)
    else:
        com = (coords * masses.reshape(-1, 1)).sum(axis=0) / total_mass

    coords_shifted = coords - com
    return coords_shifted


def _write_conf_xyz(
    mol: Chem.Mol,
    conf_id: int,
    energy_kcal: float,
    path: Path,
    move_to_com: bool,
    mol_id: str,
) -> None:
    """1 コンフォマーを XYZ として書き出す。"""
    path.parent.mkdir(parents=True, exist_ok=True)

    if move_to_com:
        coords = _move_positions_to_com(mol, conf_id)
    else:
        conf = mol.GetConformer(conf_id)
        n = mol.GetNumAtoms()
        coords = np.zeros((n, 3), dtype=float)
        for i in range(n):
            pos = conf.GetAtomPosition(i)
            coords[i] = [pos.x, pos.y, pos.z]

    lines: List[str] = []
    n_atoms = mol.GetNumAtoms()
    lines.append(str(n_atoms))
    lines.append(f"ID={mol_id} conf={conf_id} E={energy_kcal:.6f} kcal/mol")

    for i in range(n_atoms):
        atom = mol.GetAtomWithIdx(i)
        sym = atom.GetSymbol()
        x, y, z = coords[i]
        lines.append(f"{sym:2s}  {x: .8f}  {y: .8f}  {z: .8f}")

    path.write_text("\n".join(lines), encoding="utf-8")


# =========================
# 公開関数
# =========================

def run_mm_for_job(
    mol_id: str,
    smiles: str,
    embed_cfg: EmbedConfig,
    cfg: MMConfig,
    out_cfg: OutputConfig,
) -> MMResult:
    """
    1 ID について MM 構造を生成・最適化し、

      output/<ID>/MM       : すべての MM 最適化済み構造
      output/<ID>/post_MM  : エネルギー & RMSD 重複除去後の上位構造（UMA に渡す）
      output/<ID>/logs/mm_summary.log : 構造数のサマリ

    を作成する。
    """
    base_dir = out_cfg.base_dir
    id_dir = id_root_dir(base_dir, mol_id)
    mm_out_dir = mm_dir(id_dir)
    post_out_dir = post_mm_dir(id_dir)
    log_dir = logs_dir(id_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    # 1) SMILES から初期コンフォマー集合を作る
    mol, conf_ids = _prepare_initial_confs_for_smiles(smiles, embed_cfg=embed_cfg)

    # 2) MM 最適化
    energies = _mm_optimize_all_confs(mol, conf_ids, cfg)

    # 3) エネルギー昇順に並べる
    items = sorted(energies.items(), key=lambda kv: kv[1])  # (conf_id, energy) 昇順

    # 4) 重複除去＋選抜（post_MM 用）
    kept_cids = _select_unique_confs(mol, energies, items, cfg)

    # 5) XYZ 出力
    mm_out_dir.mkdir(parents=True, exist_ok=True)
    post_out_dir.mkdir(parents=True, exist_ok=True)

    full_xyz_by_cid: Dict[int, str] = {}
    all_xyz_files: List[str] = []

    # 5-1) すべて MM フォルダに書き出し
    for rank, (cid, e) in enumerate(items, start=1):
        fname = f"conf_{rank:04d}.xyz"
        fpath = Path(mm_out_dir) / fname
        _write_conf_xyz(
            mol,
            conf_id=cid,
            energy_kcal=e,
            path=fpath,
            move_to_com=out_cfg.move_to_com_origin,
            mol_id=mol_id,
        )
        full_xyz_by_cid[cid] = str(fpath)
        all_xyz_files.append(str(fpath))

    # 5-2) 選抜されたものだけ post_MM にコピー（新しい連番）
    post_xyz_files: List[str] = []
    for pr, cid in enumerate(kept_cids, start=1):
        src = Path(full_xyz_by_cid[cid])
        dst = Path(post_out_dir) / f"conf_{pr:04d}.xyz"
        shutil.copy2(src, dst)
        post_xyz_files.append(str(dst))

    # 6) mm_summary.log を ID ごとの logs/ に出力
    summary_path = Path(log_dir) / "mm_summary.log"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write(f"ID: {mol_id}\n")
        f.write(f"generated_confs: {len(conf_ids)}\n")
        f.write(f"mm_optimized:   {len(energies)}\n")
        f.write(f"unique_selected_for_UMA: {len(kept_cids)}\n")
        f.write(f"equal_tol(kcal/mol): {getattr(cfg, 'equal_tol', 0.003)}\n")
        f.write(f"energy_window(kcal/mol): {getattr(cfg, 'energy_window', 10.0)}\n")
        f.write(f"rmsd_min(Å): {getattr(cfg, 'rmsd_min', 0.8)}\n")
        f.write(f"max_keep: {getattr(cfg, 'max_keep', 10)}\n")

    return MMResult(
        n_generated=len(conf_ids),
        n_optimized=len(energies),
        xyz_files=post_xyz_files if post_xyz_files else all_xyz_files,
    )









