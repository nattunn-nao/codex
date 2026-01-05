from __future__ import annotations

from pathlib import Path
from typing import List, Dict, Optional


def _safe_float(token: str) -> Optional[float]:
    try:
        return float(token)
    except (TypeError, ValueError):
        return None


def _is_all_int_tokens(tokens: List[str]) -> bool:
    if not tokens:
        return False
    return all(t.lstrip("+-").isdigit() for t in tokens)


def collect_mo_energies_from_molecular_orbitals(lines: List[str]) -> List[float]:
    """
    最後の 'MOLECULAR ORBITALS' セクションから
    MO エネルギー (eV とみなす) を 1 次元リストで取得。
    """
    energies_ev: List[float] = []
    in_section = False

    for i, line in enumerate(lines):
        if "MOLECULAR ORBITALS" in line:
            in_section = True
            energies_ev = []  # 最後のセクションだけ
            continue

        if not in_section:
            continue

        stripped = line.strip()
        if not stripped:
            continue

        if any(ch.isalpha() for ch in stripped):
            continue

        tokens = stripped.split()

        # MO 番号行
        if _is_all_int_tokens(tokens):
            if i + 1 < len(lines):
                energy_line = lines[i + 1].strip()
                if energy_line and not any(ch.isalpha() for ch in energy_line):
                    for t in energy_line.split():
                        v = _safe_float(t)
                        if v is not None:
                            energies_ev.append(v)
            continue

    return energies_ev


def parse_occupied_mo_count(lines: List[str]) -> Optional[int]:
    """
    'NUMBER OF OCCUPIED ORBITALS' から占有MO数を取得。
    (ALPHA) があればそれを優先。
    """
    alpha_occ: Optional[float] = None
    simple_occ: Optional[int] = None

    for line in lines:
        if "NUMBER OF OCCUPIED ORBITALS" not in line:
            continue

        if "(ALPHA" in line or "(ALPHA)" in line:
            tokens = line.replace("=", " ").split()
            nums = [_safe_float(t) for t in tokens if _safe_float(t) is not None]
            if nums:
                alpha_occ = nums[-1]
        else:
            tokens = line.replace("=", " ").split()
            nums_i = [_safe_float(t) for t in tokens if _safe_float(t) is not None]
            if nums_i:
                simple_occ = int(nums_i[-1])

    if alpha_occ is not None:
        return int(alpha_occ)
    if simple_occ is not None:
        return simple_occ
    return None


def extract_frontier_mos_from_lines(lines: List[str]) -> Dict[str, Optional[float]]:
    """
    行リストから HOMO / next_HOMO / LUMO / next_LUMO (eV) を抽出。

    - HOMO            : 最高被占軌道
    - next_HOMO       : HOMO-1
    - LUMO            : 最低空軌道
    - next_LUMO       : LUMO+1
    """
    energies_ev = collect_mo_energies_from_molecular_orbitals(lines)
    if not energies_ev:
        return {
            "HOMO_eV": None,
            "next_HOMO_eV": None,
            "LUMO_eV": None,
            "next_LUMO_eV": None,
        }

    n_occ = parse_occupied_mo_count(lines)
    if n_occ is None or n_occ <= 0 or n_occ >= len(energies_ev):
        return {
            "HOMO_eV": None,
            "next_HOMO_eV": None,
            "LUMO_eV": None,
            "next_LUMO_eV": None,
        }

    homo_idx = n_occ - 1
    lumo_idx = n_occ

    homo = energies_ev[homo_idx]
    next_homo = energies_ev[homo_idx - 1] if homo_idx - 1 >= 0 else None  # HOMO-1
    lumo = energies_ev[lumo_idx]
    next_lumo = energies_ev[lumo_idx + 1] if lumo_idx + 1 < len(energies_ev) else None  # LUMO+1

    return {
        "HOMO_eV": homo,
        "next_HOMO_eV": next_homo,
        "LUMO_eV": lumo,
        "next_LUMO_eV": next_lumo,
    }


def extract_frontier_mos_from_file(path: Path) -> Dict[str, Optional[float]]:
    """
    1つの .out ファイルから frontier MO を抽出するラッパー。
    """
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except Exception as e:
        print(f"[WARN] ファイルを読み込めませんでした: {path} ({e})")
        return {
            "HOMO_eV": None,
            "next_HOMO_eV": None,
            "LUMO_eV": None,
            "next_LUMO_eV": None,
        }

    lines = text.splitlines()
    return extract_frontier_mos_from_lines(lines)
