# conformer_app/qm_layer/parser.py
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def _read_lines(path: Path) -> List[str]:
    return path.read_text(encoding="utf-8", errors="ignore").splitlines()


def _extract_final_energy(lines: List[str]) -> Optional[float]:
    """
    GAMESS の .out から最終エネルギー (Hartree) を推定する。

    優先:
      1) NSERCH: ... E= ... 行の最後の出現
      2) "FINAL ENERGY" / "TOTAL ENERGY" にマッチする行
    """
    last_e = None

    # 1) NSERCH: 行
    pat_nserch = re.compile(r"NSERCH:\s*\d+.*E=\s*([\-0-9\.]+)")
    for line in lines:
        m = pat_nserch.search(line)
        if m:
            try:
                last_e = float(m.group(1))
            except ValueError:
                pass

    if last_e is not None:
        return last_e

    # 2) "FINAL ENERGY" / "TOTAL ENERGY"
    pat_final = re.compile(r"(FINAL|TOTAL)\s+ENERGY\s*=\s*([\-0-9\.]+)")
    for line in lines:
        m = pat_final.search(line)
        if m:
            try:
                return float(m.group(2))
            except ValueError:
                continue

    return None


def _find_equilibrium_coords_header(lines: List[str]) -> Optional[int]:
    """
    'EQUILIBRIUM GEOMETRY LOCATED' の「最後の出現」以降で、
    'COORDINATES OF ALL ATOMS ARE' ブロックのヘッダ行を探す。
    """
    eq_indices: List[int] = [
        i for i, line in enumerate(lines)
        if "EQUILIBRIUM GEOMETRY LOCATED" in line
    ]
    if not eq_indices:
        return None

    start = eq_indices[-1]  # 最後の平衡構造を採用
    for j in range(start, len(lines)):
        if "COORDINATES OF ALL ATOMS ARE" in lines[j]:
            return j

    return None


def _find_last_coords_header(lines: List[str]) -> Optional[int]:
    """
    ファイル全体から最後の 'COORDINATES OF ALL ATOMS ARE' を探す。
    """
    header_idx: List[int] = []
    for i, line in enumerate(lines):
        if "COORDINATES OF ALL ATOMS ARE" in line:
            header_idx.append(i)
    if not header_idx:
        return None
    return header_idx[-1]


def _parse_coords_block_from_header(
    lines: List[str],
    header_idx: int,
) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """
    'COORDINATES OF ALL ATOMS ARE' の行インデックスから、
    直後の座標ブロックをパースする。

    対応フォーマットの例：

    (1) あなたのケース：
        COORDINATES OF ALL ATOMS ARE (ANGS)
          ATOM   CHARGE       X              Y              Z
        ------------------------------------------------------------
        N           7.0   1.5339 ...   0.4496 ...   0.2925 ...

        → tokens = [ "N", "7.0", x, y, z ]
          symbol = "N"

    (2) インデックス付き：
        COORDINATES OF ALL ATOMS ARE (ANGS)
          ATOM   ATOMIC                      COORDINATES (ANGS)
                 CHARGE         X         Y         Z
           1 C    6.0      0.000000   0.000000   0.0000z00

        → tokens = [ "1", "C", "6.0", x, y, z ]
          symbol = "C"
    """
    # ヘッダ数行をスキップするために、少し下から探索する
    i = header_idx + 1

    # まず「====/----」などの区切り線を過ぎたあたりまで進める
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if set(line.replace("-", "")) == set(""):
            # 全部 '-'
            i += 1
            continue
        # ここから下を座標ブロック候補とする
        break

    symbols: List[str] = []
    coords: List[Tuple[float, float, float]] = []

    started = False

    while i < len(lines):
        line = lines[i].strip()
        if not line:
            # 一度でも座標を読んでいれば、空行でブロック終了とみなす
            if started:
                break
            i += 1
            continue

        # 区切り線っぽい行はスキップ
        if set(line.replace("-", "")) == set(""):
            if started:
                break
            i += 1
            continue

        parts = line.split()
        if len(parts) < 4:
            # まだ座標読み込みが始まってなければ無視して進む
            if not started:
                i += 1
                continue
            # 既に座標を読み始めているなら、ここでブロック終了とみなす
            break

        # 最後の3トークンを座標として解釈できるか試す
        try:
            x = float(parts[-3])
            y = float(parts[-2])
            z = float(parts[-1])
        except ValueError:
            if not started:
                i += 1
                continue
            break

        # シンボル推定:
        #   先頭トークンが整数 → [idx, sym, ...]
        #   そうでなければ先頭トークンをシンボルとみなす
        sym_token = parts[0]
        try:
            int(sym_token)
            # ここに来たら index + symbol 形式を想定
            if len(parts) >= 2:
                symbol = parts[1]
            else:
                symbol = sym_token
        except ValueError:
            symbol = sym_token

        symbols.append(symbol)
        coords.append((x, y, z))
        started = True
        i += 1

    if not symbols:
        raise RuntimeError("座標ブロックのパースに失敗しました（原子行が見つかりません）。")

    return symbols, coords


def _extract_final_geometry(lines: List[str]) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """
    GAMESS の .out から「平衡構造の最終座標」を抜き出す。

    優先順位:
      1) 'EQUILIBRIUM GEOMETRY LOCATED' 以降で最初に現れる
         'COORDINATES OF ALL ATOMS ARE' ブロック
      2) 1) が見つからない場合は、ファイル全体で最後の
         'COORDINATES OF ALL ATOMS ARE' ブロック
    """
    # 1) EQUILIBRIUM GEOMETRY LOCATED → そこから下で座標ブロック
    header_idx = _find_equilibrium_coords_header(lines)
    if header_idx is None:
        # 2) フォールバック: 最後の COORDINATES ブロック
        header_idx = _find_last_coords_header(lines)

    if header_idx is None:
        raise RuntimeError(
            "GAMESS 出力中に 'EQUILIBRIUM GEOMETRY LOCATED' も "
            "'COORDINATES OF ALL ATOMS ARE' も見つかりません。"
        )

    return _parse_coords_block_from_header(lines, header_idx)


def _write_xyz(path: Path, symbols: List[str], coords: List[Tuple[float, float, float]], comment: str = "") -> None:
    """非常に簡単な XYZ writer"""
    n = len(symbols)
    lines: List[str] = []
    lines.append(str(n))
    lines.append(comment)
    for s, (x, y, z) in zip(symbols, coords):
        lines.append(f"{s:2s}  {x: .8f}  {y: .8f}  {z: .8f}")
    path.write_text("\n".join(lines), encoding="utf-8")


def extract_final_structure_and_energy(
    job_dir: str,
    slot: str,
    out_basename: Optional[str] = None,
) -> Dict[str, object]:
    """
    1 つの gmsX フォルダについて、
      - gmsX.out をパースして最終エネルギーと「平衡構造の最終座標」を取得
      - job_dir/<out_basename or 'final.xyz'> に XYZ を出力
    """
    job_dir_path = Path(job_dir)
    out_file = job_dir_path / f"{slot}.out"
    if out_basename is None:
        xyz_path = job_dir_path / "final.xyz"
    else:
        xyz_path = job_dir_path / out_basename

    if not out_file.exists():
        return {
            "slot": slot,
            "job_dir": str(job_dir_path),
            "out_path": str(out_file),
            "status": "error",
            "energy_hartree": None,
            "final_xyz": "",
            "message": f"出力ファイルが見つかりません: {out_file}",
        }

    lines = _read_lines(out_file)

    # エネルギー
    try:
        energy = _extract_final_energy(lines)
    except Exception as e:
        return {
            "slot": slot,
            "job_dir": str(job_dir_path),
            "out_path": str(out_file),
            "status": "error",
            "energy_hartree": None,
            "final_xyz": "",
            "message": f"エネルギー取得に失敗: {e}",
        }

    # 平衡構造の座標
    try:
        symbols, coords = _extract_final_geometry(lines)
        comment = f"{slot}  E={energy if energy is not None else 'NA'} Hartree"
        _write_xyz(xyz_path, symbols, coords, comment=comment)
    except Exception as e:
        return {
            "slot": slot,
            "job_dir": str(job_dir_path),
            "out_path": str(out_file),
            "status": "error",
            "energy_hartree": energy,
            "final_xyz": "",
            "message": f"最終構造の取得または XYZ 書き出しに失敗: {e}",
        }

    return {
        "slot": slot,
        "job_dir": str(job_dir_path),
        "out_path": str(out_file),
        "status": "ok",
        "energy_hartree": energy,
        "final_xyz": str(xyz_path),
        "message": "",
    }
