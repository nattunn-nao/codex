from __future__ import annotations

"""
extract.py (geom_qm_extract_logic)

指定したフォルダ内の

  - GAMESS .out ファイル
  - *_final.xyz ファイル

を自動検出し、それぞれに対して

  - .out  : 双極子モーメント + frontier MO (HOMO/LUMO) を抽出
  - .xyz : SASA / vdW surface / vdW volume / Rg / Ovality を計算

した結果をマージし、Geom_QM_summary.xlsx にまとめて出力するロジックモジュール。

対応関係の例：
  lmp.out       → id = "lmp"
  lmp_final.xyz → id = "lmp"

想定配置（例）:
  - conformer_app/
      Geom_module/
          geom_props_logic.py      : analyze_xyz_file()
          sasa_calc_logic.py       : PROBE_RADIUS, GRID_SPACING
      QM_module/
          gamess_dipole.py         : parse_dipole_from_lines(lines)
          gamess_mo.py             : extract_frontier_mos_from_lines(lines)
      extract/
          extract.py               : いまこのファイル

メインエントリ（外部から呼び出す想定）:
    run_geom_qm_extract(root: Path) -> Path | None

CLI としても利用可能:
    python -m conformer_app.extract.extract <対象フォルダパス>
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Any

import pandas as pd

from .geom_props_logic import analyze_xyz_file
from .sasa_calc_logic import PROBE_RADIUS, GRID_SPACING
from .gamess_dipole import parse_dipole_from_lines
from .gamess_mo import extract_frontier_mos_from_lines


# ============================================================
# ユーティリティ
# ============================================================

def make_id_from_out(path: Path) -> str:
    """
    .out 用のID。
    例: sample.out → "sample"
    """
    return path.stem


def make_id_from_xyz(path: Path) -> str:
    """
    *_final.xyz 用のID。

    例:
      sample_final.xyz → "sample"
      foo.xyz          → "foo"   （_final で終わっていなければそのまま stem）
    """
    stem = path.stem
    if stem.endswith("_final"):
        return stem[:-6]  # remove "_final"
    return stem


# ============================================================
# QM (.out) 側
# ============================================================

def process_gamess_out_file(path: Path) -> Dict[str, Optional[float]]:
    """
    1 つの GAMESS .out から双極子 + frontier MO を抽出して dict を返す。

    戻り値キー:
      - id
      - out_filename
      - dipole_x, dipole_y, dipole_z, dipole_moment
      - HOMO_eV, next_HOMO_eV, LUMO_eV, next_LUMO_eV
    """
    rec: Dict[str, Optional[float]] = {
        "id": make_id_from_out(path),
        "out_filename": path.name,
        "dipole_x": None,
        "dipole_y": None,
        "dipole_z": None,
        "dipole_moment": None,
        "HOMO_eV": None,
        "next_HOMO_eV": None,
        "LUMO_eV": None,
        "next_LUMO_eV": None,
    }

    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except Exception as e:
        print(f"[WARN] .out ファイルを読み込めませんでした: {path} ({e})")
        return rec

    lines = text.splitlines()

    dipole = parse_dipole_from_lines(lines)
    mos = extract_frontier_mos_from_lines(lines)

    # dipole 部分
    rec["dipole_x"] = dipole.get("dipole_x")
    rec["dipole_y"] = dipole.get("dipole_y")
    rec["dipole_z"] = dipole.get("dipole_z")
    rec["dipole_moment"] = dipole.get("dipole_moment")

    # frontier MO 部分
    rec["HOMO_eV"] = mos.get("HOMO_eV")
    rec["next_HOMO_eV"] = mos.get("next_HOMO_eV")
    rec["LUMO_eV"] = mos.get("LUMO_eV")
    rec["next_LUMO_eV"] = mos.get("next_LUMO_eV")

    return rec


def collect_qm_records(dir_path: Path) -> Dict[str, Dict[str, Any]]:
    """
    dir_path 内の *.out をすべて処理し、
    id → レコード dict の辞書を返す。
    """
    qm_files = sorted(dir_path.glob("*.out"))
    if not qm_files:
        print("[QM] 対象となる .out ファイルがありません。")
        return {}

    print(f"[QM] {len(qm_files)} 個の .out ファイルを処理します。")
    id_to_rec: Dict[str, Dict[str, Any]] = {}

    for p in qm_files:
        print(f"  [QM] 処理中: {p.name}")
        rec = process_gamess_out_file(p)
        id_to_rec[rec["id"]] = rec

    return id_to_rec


# ============================================================
# Geom (.xyz) 側
# ============================================================

def process_xyz_file(path: Path) -> Dict[str, Any]:
    """
    1 つの *_final.xyz から幾何プロパティを計算して dict を返す。

    戻り値キー:
      - id
      - xyz_filename
      - n_atoms, n_polar_atoms
      - probe_radius_A, grid_spacing_A
      - SASA_area_A2, SASA_polar_area_A2
      - vdW_area_A2, vdW_volume_A3
      - Rg_A, ovality
    """
    res = analyze_xyz_file(
        path,
        write_html=False,   # HTML は出さない
        write_excel=False,  # per-XYZ Excel もここでは作らない
    )

    polar_mask = res["polar_mask"]

    rec: Dict[str, Any] = {
        "id": make_id_from_xyz(path),
        "xyz_filename": Path(res["xyz_path"]).name,
        "n_atoms": len(res["symbols"]),
        "n_polar_atoms": int(polar_mask.sum()),
        "probe_radius_A": PROBE_RADIUS,
        "grid_spacing_A": GRID_SPACING,
        "SASA_area_A2": res["sasa_area"],
        "SASA_polar_area_A2": res["sasa_polar_area"],
        "vdW_area_A2": res["vdw_area"],
        "vdW_volume_A3": res["vdw_volume"],
        "Rg_A": res["Rg"],
        "ovality": res["ovality"],
    }
    return rec


def collect_geom_records(dir_path: Path) -> Dict[str, Dict[str, Any]]:
    """
    dir_path 内の *_final.xyz をすべて処理し、
    id → レコード dict の辞書を返す。
    """
    xyz_files = sorted(dir_path.glob("*_final.xyz"))
    if not xyz_files:
        print("[Geom] 対象となる *_final.xyz ファイルがありません。")
        return {}

    print(f"[Geom] {len(xyz_files)} 個の *_final.xyz ファイルを処理します。")
    id_to_rec: Dict[str, Dict[str, Any]] = {}

    for p in xyz_files:
        print(f"  [Geom] 処理中: {p.name}")
        rec = process_xyz_file(p)
        id_to_rec[rec["id"]] = rec

        # コンソールには軽くサマリだけ
        print(f"    原子数: {rec['n_atoms']}")
        print(f"    極性原子数 (heavy + H): {rec['n_polar_atoms']}")
        print(f"    SASA (Å^2): {rec['SASA_area_A2']:.3f}")
        print(f"    SASA (polar, Å^2): {rec['SASA_polar_area_A2']:.3f}")
        print(f"    vdW surface area (Å^2): {rec['vdW_area_A2']:.3f}")
        print(f"    vdW volume (Å^3): {rec['vdW_volume_A3']:.3f}")
        print(f"    Radius of gyration Rg (Å): {rec['Rg_A']:.3f}")
        print(f"    Ovality (from vdW S, V): {rec['ovality']:.3f}")

    return id_to_rec


# ============================================================
# QM + Geom のマージ & Excel 出力
# ============================================================

def merge_qm_geom(
    qm_records: Dict[str, Dict[str, Any]],
    geom_records: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    QM 側と Geom 側のレコードを id でマージし、
    行のリストを返す。
    """
    all_ids = sorted(set(qm_records.keys()) | set(geom_records.keys()))
    rows: List[Dict[str, Any]] = []

    for _id in all_ids:
        qm = qm_records.get(_id, {})
        gm = geom_records.get(_id, {})

        row: Dict[str, Any] = {"id": _id}

        # ファイル名
        row["out_filename"] = qm.get("out_filename")
        row["xyz_filename"] = gm.get("xyz_filename")

        # QM 部分
        row["dipole_x"] = qm.get("dipole_x")
        row["dipole_y"] = qm.get("dipole_y")
        row["dipole_z"] = qm.get("dipole_z")
        row["dipole_moment"] = qm.get("dipole_moment")

        row["HOMO_eV"] = qm.get("HOMO_eV")
        row["next_HOMO_eV"] = qm.get("next_HOMO_eV")
        row["LUMO_eV"] = qm.get("LUMO_eV")
        row["next_LUMO_eV"] = qm.get("next_LUMO_eV")

        # Geom 部分
        row["n_atoms"] = gm.get("n_atoms")
        row["n_polar_atoms"] = gm.get("n_polar_atoms")
        row["probe_radius_A"] = gm.get("probe_radius_A")
        row["grid_spacing_A"] = gm.get("grid_spacing_A")
        row["SASA_area_A2"] = gm.get("SASA_area_A2")
        row["SASA_polar_area_A2"] = gm.get("SASA_polar_area_A2")
        row["vdW_area_A2"] = gm.get("vdW_area_A2")
        row["vdW_volume_A3"] = gm.get("vdW_volume_A3")
        row["Rg_A"] = gm.get("Rg_A")
        row["ovality"] = gm.get("ovality")

        rows.append(row)

    return rows


def write_combined_excel(rows: List[Dict[str, Any]], out_path: Path) -> None:
    df = pd.DataFrame(rows)
    df.to_excel(out_path, index=False)
    print(f"[ALL] 統合 Excel 出力: {out_path}")


# ============================================================
# 外部呼び出し用エントリポイント
# ============================================================

def run_geom_qm_extract(root: Path) -> Optional[Path]:
    """
    外部コード（QM ポスト処理など）から呼び出す想定の関数。

    Parameters
    ----------
    root : Path
        .out と *_final.xyz が格納されているディレクトリ。
        例: result/ID フォルダ

    Returns
    -------
    out_path : Path | None
        Geom_QM_summary.xlsx のパス。
        解析対象が見つからなかった場合は None。
    """
    root = root.resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"指定されたパスはディレクトリではありません: {root}")

    print("========================================")
    print(f"Geom+QM 統合解析 対象フォルダ: {root}")
    print("========================================")

    # QM 側
    qm_records = collect_qm_records(root)

    # Geom 側
    geom_records = collect_geom_records(root)

    if not qm_records and not geom_records:
        print("[ALL] .out も *_final.xyz も見つからなかったので、何も出力しません。")
        return None

    # マージ & Excel 出力
    rows = merge_qm_geom(qm_records, geom_records)
    out_xlsx = root / "Geom_QM_summary.xlsx"
    write_combined_excel(rows, out_xlsx)

    print("========================================")
    print("Geom+QM 統合解析 完了")
    print("========================================")

    return out_xlsx


# ============================================================
# CLI 本体
# ============================================================

def main(argv: List[str] | None = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    if not argv:
        print("使い方: python -m conformer_app.extract.extract <対象フォルダパス>")
        print("  例:   python -m conformer_app.extract.extract D:\\path\\to\\result\\ID1")
        sys.exit(1)

    root = Path(argv[0])
    try:
        out_path = run_geom_qm_extract(root)
        if out_path is None:
            sys.exit(0)
    except Exception as e:
        print(f"[ERROR] 処理中にエラーが発生しました: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
