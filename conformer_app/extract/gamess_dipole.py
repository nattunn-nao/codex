from __future__ import annotations

from pathlib import Path
from typing import List, Dict, Optional


def _safe_float(token: str) -> Optional[float]:
    try:
        return float(token)
    except (TypeError, ValueError):
        return None


def _parse_system_electrostatic_block(
    lines: List[str],
    start_idx: int,
    max_len: int = 80,
) -> Dict[str, Optional[float]]:
    """
    SYSTEM ELECTROSTATIC MOMENTS ブロックから双極子を抽出。
    可能なら PCM 表面電荷込みの TOTAL DIPOLE を優先。
    """
    result: Dict[str, Optional[float]] = {
        "dipole_x": None,
        "dipole_y": None,
        "dipole_z": None,
        "dipole_moment": None,  # total (Debye)
    }

    end_idx = min(len(lines), start_idx + max_len)
    seen_pcm_total = False
    last_header_idx: Optional[int] = None

    for i in range(start_idx, end_idx):
        line = lines[i]

        if "INCLUDING PCM SURFACE CHARGES" in line and "TOTAL DIPOLE IS" in line:
            seen_pcm_total = True
            continue

        if (
            "DX" in line
            and "DY" in line
            and "DZ" in line
            and ("/D/" in line or "DEBYE" in line)
        ):
            # PCM TOTAL 優先
            if seen_pcm_total:
                for j in range(i + 1, end_idx):
                    row = lines[j].strip()
                    if not row:
                        continue
                    vals = [_safe_float(t) for t in row.split() if _safe_float(t) is not None]
                    if len(vals) >= 4:
                        result["dipole_x"], result["dipole_y"], result["dipole_z"], result["dipole_moment"] = vals[:4]
                        return result
                    break
            # fallback 用
            last_header_idx = i

    # PCM TOTAL が取れなかった場合
    if last_header_idx is not None:
        for j in range(last_header_idx + 1, end_idx):
            row = lines[j].strip()
            if not row:
                continue
            vals = [_safe_float(t) for t in row.split() if _safe_float(t) is not None]
            if len(vals) >= 4:
                result["dipole_x"], result["dipole_y"], result["dipole_z"], result["dipole_moment"] = vals[:4]
                return result

    return result


def parse_dipole_from_lines(lines: List[str]) -> Dict[str, Optional[float]]:
    """
    GAMESS 出力（行リスト）から双極子モーメントを抽出する。
    戻り値: {"dipole_x", "dipole_y", "dipole_z", "dipole_moment"}
    """
    default: Dict[str, Optional[float]] = {
        "dipole_x": None,
        "dipole_y": None,
        "dipole_z": None,
        "dipole_moment": None,
    }

    # 1) SYSTEM ELECTROSTATIC MOMENTS の最後のブロック
    last_system_idx: Optional[int] = None
    for i, line in enumerate(lines):
        if "SYSTEM ELECTROSTATIC MOMENTS" in line:
            last_system_idx = i

    if last_system_idx is not None:
        d = _parse_system_electrostatic_block(lines, last_system_idx)
        if any(v is not None for v in d.values()):
            return d

    # 2) フォールバック: DIPOLE MOMENT (DEBYE)
    dipole = default.copy()
    in_block = False

    for line in lines:
        if "DIPOLE MOMENT" in line and "DEBYE" in line:
            in_block = True
            dipole = default.copy()
            continue

        if not in_block:
            continue

        row = line.strip()
        if not row:
            continue

        tokens = row.split()

        # X 0.0 / TOTAL 1.0 形式
        if tokens[0] in ("X", "Y", "Z", "TOTAL") and len(tokens) >= 2:
            label = tokens[0]
            value = _safe_float(tokens[1])
            if label == "X":
                dipole["dipole_x"] = value
            elif label == "Y":
                dipole["dipole_y"] = value
            elif label == "Z":
                dipole["dipole_z"] = value
            elif label == "TOTAL":
                dipole["dipole_moment"] = value
                in_block = False
            continue

        # X 0.0 Y 0.1 Z 0.2 TOTAL 0.3 形式
        for label in ("X", "Y", "Z", "TOTAL"):
            if label in tokens:
                idx = tokens.index(label)
                if idx + 1 < len(tokens):
                    value = _safe_float(tokens[idx + 1])
                    if label == "X":
                        dipole["dipole_x"] = value
                    elif label == "Y":
                        dipole["dipole_y"] = value
                    elif label == "Z":
                        dipole["dipole_z"] = value
                    elif label == "TOTAL":
                        dipole["dipole_moment"] = value
                        in_block = False

    return dipole


def extract_dipole_from_file(path: Path) -> Dict[str, Optional[float]]:
    """
    1つの .out ファイルから双極子を抽出するラッパー。
    """
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except Exception as e:
        print(f"[WARN] ファイルを読み込めませんでした: {path} ({e})")
        return {
            "dipole_x": None,
            "dipole_y": None,
            "dipole_z": None,
            "dipole_moment": None,
        }

    lines = text.splitlines()
    return parse_dipole_from_lines(lines)
