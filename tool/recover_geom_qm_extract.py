from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Any

import pandas as pd

"""
python .\tool\recover_geom_qm_extract.py .\result\ID\ID_best.xyz
"""

HERE = Path(__file__).resolve()
PROJECT_ROOT = HERE.parent.parent  # re_production/
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# ここから conformer_app をインポート
from conformer_app.extract.geom_props_logic import analyze_xyz_file
from conformer_app.extract.sasa_calc_logic import PROBE_RADIUS, GRID_SPACING
from conformer_app.extract.gamess_dipole import parse_dipole_from_lines
from conformer_app.extract.gamess_mo import extract_frontier_mos_from_lines


# ============================================================
# ユーティリティ
# ============================================================


def _infer_id_from_best_xyz(xyz_path: Path) -> str:
    """
    4_best.xyz → ID = "4"
    111-1_best.xyz → ID = "111-1"
    などを想定。
    _best が無いときは stem 全体を ID とみなし、
    それも変な場合は親フォルダ名を ID とする。
    """
    stem = xyz_path.stem  # 例: "4_best"
    if stem.endswith("_best"):
        return stem[:-5]
    # 念のため空文字は避ける
    if stem:
        return stem
    return xyz_path.parent.name


def _find_qm_out_for_best(id_str: str) -> Path:
    """
    output/<ID>/QM/gms1/gms1.out を第一候補とし、
    なければ output/<ID>/QM/gms*/gms*.out のうち一つを返す。
    見つからなければ RuntimeError。
    """
    output_root = PROJECT_ROOT / "output"
    qm_root = output_root / id_str / "QM"

    if not qm_root.is_dir():
        raise RuntimeError(f"QM ディレクトリが見つかりません: {qm_root}")

    # 1. gms1/gms1.out を優先
    gms1_dir = qm_root / "gms1"
    gms1_out = gms1_dir / "gms1.out"
    if gms1_out.is_file():
        return gms1_out

    # 2. gms*/gms*.out のうち、最初に見つかったもの
    for sub in sorted(qm_root.glob("gms*")):
        if not sub.is_dir():
            continue
        cand = sub / f"{sub.name}.out"  # gms2/gms2.out など
        if cand.is_file():
            return cand

    raise RuntimeError(f"QM .out ファイルが見つかりません: {qm_root}")


def _process_qm_out(out_path: Path) -> Dict[str, Any]:
    """
    GAMESS .out を読み込み、双極子 + frontier MO を抽出して dict で返す。
    """
    text = out_path.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()

    dipole = parse_dipole_from_lines(lines)
    mos = extract_frontier_mos_from_lines(lines)

    rec: Dict[str, Any] = {
        "out_filename": out_path.name,
        "dipole_x": dipole.get("dipole_x"),
        "dipole_y": dipole.get("dipole_y"),
        "dipole_z": dipole.get("dipole_z"),
        "dipole_moment": dipole.get("dipole_moment"),
        "HOMO_eV": mos.get("HOMO_eV"),
        "next_HOMO_eV": mos.get("next_HOMO_eV"),
        "LUMO_eV": mos.get("LUMO_eV"),
        "next_LUMO_eV": mos.get("next_LUMO_eV"),
    }
    return rec


# ============================================================
# メイン処理
# ============================================================


def recover_for_best_xyz(best_xyz_path: Path) -> Path:
    """
    4_best.xyz のような XYZ を起点に

      1) XYZ から Geom (SASA / vdW / Rg / Ovality) を再計算
      2) output/<ID>/QM/gms1/gms1.out などから QM 情報を抽出
      3) result/<ID>/Geom_QM_summary.xlsx を新規作成 or 上書き

    を行い、その Excel パスを返す。
    """
    best_xyz_path = best_xyz_path.resolve()
    if not best_xyz_path.is_file():
        raise FileNotFoundError(f"best XYZ が見つかりません: {best_xyz_path}")

    # ID の推定
    mol_id = _infer_id_from_best_xyz(best_xyz_path)

    # 1) XYZ 解析（HTML/Excel はここでは出さない）
    geom_res = analyze_xyz_file(
        best_xyz_path,
        write_html=False,
        write_excel=False,
    )

    # Geom 部分の値を取り出し
    symbols = geom_res["symbols"]
    polar_mask = geom_res["polar_mask"]
    rec_geom: Dict[str, Any] = {
        "id": mol_id,
        "xyz_filename": Path(geom_res["xyz_path"]).name,
        "n_atoms": len(symbols),
        "n_polar_atoms": int(polar_mask.sum()),
        "probe_radius_A": PROBE_RADIUS,
        "grid_spacing_A": GRID_SPACING,
        "SASA_area_A2": geom_res["sasa_area"],
        "SASA_polar_area_A2": geom_res["sasa_polar_area"],
        "vdW_area_A2": geom_res["vdw_area"],
        "vdW_volume_A3": geom_res["vdw_volume"],
        "Rg_A": geom_res["Rg"],
        "ovality": geom_res["ovality"],
    }

    # 2) 対応する QM .out を特定して解析
    out_path = _find_qm_out_for_best(mol_id)
    rec_qm = _process_qm_out(out_path)

    # 統合レコードを作成
    row: Dict[str, Any] = {"id": mol_id}
    row.update(rec_qm)
    row.update(rec_geom)

    df = pd.DataFrame([row])

    # 3) result/<ID>/Geom_QM_summary.xlsx に出力
    result_root = PROJECT_ROOT / "result"
    result_id_dir = result_root / mol_id
    result_id_dir.mkdir(parents=True, exist_ok=True)

    out_xlsx = result_id_dir / "Geom_QM_summary.xlsx"
    df.to_excel(out_xlsx, index=False)

    return out_xlsx


# ============================================================
# CLI エントリポイント
# ============================================================


def main(argv: list[str] | None = None) -> None:
    import sys as _sys

    if argv is None:
        argv = _sys.argv[1:]

    if not argv:
        print("使い方: python tool/recover_geom_qm.py result\\4\\4_best.xyz")
        raise SystemExit(1)

    best_xyz = Path(argv[0])

    try:
        out_xlsx = recover_for_best_xyz(best_xyz)
    except Exception as e:
        print(f"[ERROR] リカバリ中にエラーが発生しました: {e}")
        raise SystemExit(1)

    print(f"[OK] 出力ファイル: {out_xlsx}")


if __name__ == "__main__":
    main()
