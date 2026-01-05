from __future__ import annotations

import os
from pathlib import Path
from typing import Tuple, List, Dict, Optional

import numpy as np
import pandas as pd

from rdkit.Chem import rdchem, rdDetermineBonds
from rdkit.Geometry import Point3D

from .sasa_calc_logic import (
    PROBE_RADIUS,
    GRID_SPACING,
    make_grid,
    distance_field_on_grid,
    field_to_vertices_faces,
    volume_from_field,
    compute_sasa_vdw_properties,
)
from .sasa_plot_logic import plot_wireframe_colored

PT = rdchem.GetPeriodicTable()

# 極性 heavy 原子
POLAR_HEAVY = {"N", "O", "F", "P", "S", "Cl", "Br", "I"}

# 元素記号の正規化（大文字・小文字吸収）
_ELEMENT_ALIAS: dict[str, str] = {
    "SI": "Si",
    "CL": "Cl",
    "BR": "Br",
    "NA": "Na",
    "CA": "Ca",
    "MG": "Mg",
    "AL": "Al",
}


def normalize_element_symbol(sym_raw: str) -> str:
    """
    XYZ 等から読んだ元素記号を「RDKit が期待する形」に正規化する。

    - 前後の空白を除去
    - 全大文字は別名テーブルで補正 (SI -> Si, CL -> Cl, ...)
    - それ以外は 先頭大文字 + 以降小文字 (c -> C, na -> Na, PT -> Pt)
    """
    s = (sym_raw or "").strip()
    if not s:
        return s

    up = s.upper()
    if up in _ELEMENT_ALIAS:
        return _ELEMENT_ALIAS[up]

    if len(s) == 1:
        # 単一文字元素（C, N, O, F...）
        return s.upper()

    # それ以外は 先頭だけ大文字、残り小文字
    return s[0].upper() + s[1:].lower()


# ============================================================
# XYZ 読み込み
# ============================================================

def read_xyz(path: str | Path) -> Tuple[str, List[str], np.ndarray]:
    """少し頑丈な XYZ リーダー。

    - ファイル全体を読み込んで
    - 最初に int() に変換できる行を「原子数」とみなす
    """
    path = str(path)
    with open(path, "r", encoding="utf-8") as f:
        raw_lines = [l.rstrip("\n") for l in f]

    # 空行は除去
    lines = [l.strip() for l in raw_lines if l.strip()]

    if len(lines) < 2:
        raise ValueError(f"XYZ file too short: {path}")

    # 原子数行を探す
    atom_count_idx: Optional[int] = None
    n_atoms: Optional[int] = None

    for i, line in enumerate(lines):
        try:
            n_atoms = int(line)
            atom_count_idx = i
            break
        except ValueError:
            continue

    if atom_count_idx is None or n_atoms is None:
        raise ValueError(
            f"Could not find atom count in XYZ file: {path}\n"
            f"  first non-empty line = {lines[0]!r}"
        )

    # コメント行（なければファイル名）
    if atom_count_idx + 1 < len(lines):
        comment = lines[atom_count_idx + 1] or os.path.basename(path)
    else:
        comment = os.path.basename(path)

    # 原子座標行
    start = atom_count_idx + 2
    end = start + n_atoms
    if end > len(lines):
        raise ValueError(
            f"Not enough atom lines in XYZ file: {path}\n"
            f"  expected {n_atoms}, but only {len(lines) - start} lines available"
        )

    symbols: List[str] = []
    coords: List[List[float]] = []

    for line in lines[start:end]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Bad XYZ atom line in {path!r}: {line!r}")
        s_raw, x, y, z = parts[:4]
        s = normalize_element_symbol(s_raw)  # ★ ここで正規化
        symbols.append(s)
        coords.append([float(x), float(y), float(z)])

    if len(symbols) != n_atoms:
        raise ValueError(
            f"Atom count mismatch in XYZ file: {path}\n"
            f"  header says {n_atoms}, but parsed {len(symbols)} atoms"
        )

    return comment, symbols, np.array(coords, float)


# ============================================================
# RDKit Mol & 極性判定
# ============================================================

def build_mol_from_xyz(symbols: List[str], coords: np.ndarray) -> rdchem.Mol:
    """XYZ から RDKit Mol を作り、DetermineBonds で結合推定"""
    rw = rdchem.RWMol()
    conf = rdchem.Conformer(len(symbols))

    for i, (sym, (x, y, z)) in enumerate(zip(symbols, coords)):
        idx = rw.AddAtom(rdchem.Atom(sym))
        conf.SetAtomPosition(idx, Point3D(float(x), float(y), float(z)))

    rw.AddConformer(conf)
    mol = rw.GetMol()
    rdDetermineBonds.DetermineBonds(mol, useHueckel=True)
    return mol


def build_polar_mask(symbols: List[str], coords: np.ndarray) -> np.ndarray:
    """
    極性原子マスクを返す（True/False 配列）

    極性とみなす条件:
      - heavy atom が N,O,F,P,S,Cl,Br,I
      - H で、その結合相手が上記 heavy polar
    """
    mol = build_mol_from_xyz(symbols, coords)
    mask = np.zeros(len(symbols), dtype=bool)

    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        sym = atom.GetSymbol()

        if sym in POLAR_HEAVY:
            mask[i] = True
        elif sym == "H":
            if any(nb.GetSymbol() in POLAR_HEAVY for nb in atom.GetNeighbors()):
                mask[i] = True

    return mask


# ============================================================
# vdW 半径・質量・慣性半径
# ============================================================

def get_vdw_radii(symbols: List[str]) -> np.ndarray:
    """RDKit の周期表から vdW 半径を取得"""
    radii = []
    for s in symbols:
        Z = PT.GetAtomicNumber(s)
        if Z == 0:
            raise ValueError(f"Unknown element: {s}")
        r = PT.GetRvdw(Z)
        if r <= 0:
            raise ValueError(f"No vdW radius for {s}")
        radii.append(r)
    return np.array(radii, float)


def get_masses(symbols: List[str]) -> np.ndarray:
    """原子質量（原子量）を取得"""
    masses = []
    for s in symbols:
        Z = PT.GetAtomicNumber(s)
        if Z == 0:
            raise ValueError(f"Unknown element: {s}")
        m = PT.GetAtomicWeight(Z)
        masses.append(m)
    return np.array(masses, float)


def radius_of_gyration(symbols: List[str], coords: np.ndarray) -> float:
    """
    質量重み付き慣性半径 Rg (Å) を計算
    Rg^2 = sum_i m_i |r_i - r_COM|^2 / sum_i m_i
    """
    coords = np.asarray(coords, float)
    masses = get_masses(symbols)
    total_mass = masses.sum()
    if total_mass <= 0:
        masses = np.ones_like(masses)
        total_mass = masses.sum()

    com = np.sum(coords * masses[:, None], axis=0) / total_mass
    rel = coords - com
    rg2 = np.sum(masses * np.sum(rel**2, axis=1)) / total_mass
    return float(np.sqrt(rg2))


# ============================================================
# 卵形度（Ovality）
# ============================================================

def ovality_from_surface_volume(S: float, V: float) -> float:
    """
    卵形度 Ovality = S / S_min
    S_min = 4π(3V/4π)^(2/3) （同じ体積を持つ真球の表面積）

    S: 分子表面積 (Å^2)
    V: 分子体積   (Å^3)
    """
    if S <= 0 or V <= 0:
        return np.nan
    sphere_area_min = 4.0 * np.pi * (3.0 * V / (4.0 * np.pi)) ** (2.0 / 3.0)
    return float(S / sphere_area_min)


# ============================================================
# Excel 出力
# ============================================================

def save_props_to_excel(xyz_path: str,
                        n_atoms: int,
                        n_polar: int,
                        sasa_area: float,
                        sasa_polar_area: float,
                        vdw_area: float,
                        vdw_volume: float,
                        rg: float,
                        ovality: float,
                        out_xlsx: str) -> None:
    """
    分子の表面積・体積・慣性半径・卵形度を Excel に 1 行追記（存在しなければ新規作成）
    """
    row = {
        "xyz_file": os.path.basename(xyz_path),
        "n_atoms": n_atoms,
        "n_polar_atoms": n_polar,
        "probe_radius_A": PROBE_RADIUS,
        "grid_spacing_A": GRID_SPACING,
        "SASA_area_A2": sasa_area,
        "SASA_polar_area_A2": sasa_polar_area,
        "vdW_area_A2": vdw_area,
        "vdW_volume_A3": vdw_volume,
        "Rg_A": rg,
        "ovality": ovality,
    }

    if os.path.exists(out_xlsx):
        old = pd.read_excel(out_xlsx)
        df = pd.concat([old, pd.DataFrame([row])], ignore_index=True)
    else:
        df = pd.DataFrame([row])

    df.to_excel(out_xlsx, index=False)


# ============================================================
# 高レベル：XYZ 1件をまるっと解析
# ============================================================

def analyze_xyz_file(xyz_path: str | Path,
                     write_html: bool = True,
                     write_excel: bool = True) -> Dict[str, object]:
    """
    XYZ ファイル 1 つに対して、読み込み・SASA/vdW 計算・可視化・Excel 出力まで行い、
    結果を dict で返す高レベル関数。

    Parameters
    ----------
    xyz_path : str | Path
        入力 XYZ ファイルパス。
    write_html : bool, default True
        True の場合、SASA ワイヤーフレームの HTML を出力。
    write_excel : bool, default True
        True の場合、surface/volume Excel を出力。

    Returns
    -------
    result : dict
        解析結果とファイルパスなどを含む辞書。
    """
    xyz_path = str(xyz_path)
    comment, symbols, coords = read_xyz(xyz_path)

    polar_mask = build_polar_mask(symbols, coords)
    radii_vdw = get_vdw_radii(symbols)
    radii_sasa = radii_vdw + PROBE_RADIUS

    # --- グリッドと距離場（SASA 用グリッドをベースに vdW も計算） ---
    xs, ys, zs, X, Y, Z = make_grid(coords, radii_sasa, spacing=GRID_SPACING)
    field_sasa = distance_field_on_grid(X, Y, Z, coords, radii_sasa)
    field_vdw = distance_field_on_grid(X, Y, Z, coords, radii_vdw)

    # --- SASA / vdW 物性値 ---
    sasa_vdw = compute_sasa_vdw_properties(
        coords=coords,
        polar_mask=polar_mask,
        xs=xs, ys=ys, zs=zs,
        field_sasa=field_sasa,
        field_vdw=field_vdw,
    )

    # --- Rg, Ovality (vdW S, V から) ---
    rg = radius_of_gyration(symbols, coords)
    ovality = ovality_from_surface_volume(
        sasa_vdw["vdw_area"],
        sasa_vdw["vdw_volume"],
    )

    base, _ = os.path.splitext(xyz_path)
    html_path = base + "_sasa_wireframe_polar_on_sasa_with_props.html"
    xlsx_path = base + "_surface_volume.xlsx"

    # --- 可視化（SASA ワイヤーフレーム＋原子） ---
    if write_html:
        plot_wireframe_colored(
            comment=comment,
            symbols=symbols,
            coords=coords,
            verts_sasa=sasa_vdw["verts_sasa"],
            faces_sasa=sasa_vdw["faces_sasa"],
            vertex_polar=sasa_vdw["vertex_polar_sasa"],
            out_html=html_path,
        )
    else:
        html_path = None

    # --- Excel 出力 ---
    if write_excel:
        save_props_to_excel(
            xyz_path=xyz_path,
            n_atoms=len(symbols),
            n_polar=int(polar_mask.sum()),
            sasa_area=sasa_vdw["sasa_area"],
            sasa_polar_area=sasa_vdw["sasa_polar_area"],
            vdw_area=sasa_vdw["vdw_area"],
            vdw_volume=sasa_vdw["vdw_volume"],
            rg=rg,
            ovality=ovality,
            out_xlsx=xlsx_path,
        )
    else:
        xlsx_path = None

    # --- 結果をまとめて返す ---
    result: Dict[str, object] = {
        "xyz_path": xyz_path,
        "comment": comment,
        "symbols": symbols,
        "coords": coords,
        "polar_mask": polar_mask,
        "radii_vdw": radii_vdw,
        "radii_sasa": radii_sasa,
        "html_path": html_path,
        "excel_path": xlsx_path,
        "Rg": rg,
        "ovality": ovality,
    }
    result.update(sasa_vdw)
    return result
