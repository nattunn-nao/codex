# conformer_app/qm_layer/__init__.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple, Optional, Union, Callable
import re
import shutil
import subprocess

import pandas as pd

from conformer_app.qm_layer.input_builder import make_inp_for_slots
from conformer_app.uma_layer import UMAResult  # UMA 結果オブジェクト用


# ============================================================
# 型定義
# ============================================================

@dataclass
class QMResult:
    """1つの ID に対する QM 計算のまとめ結果。"""
    id: str
    n_jobs: int
    summary_path: str
    best_energy_h: Optional[float]
    best_xyz_path: str | None


# ============================================================
# ルート / template 解決ユーティリティ
# ============================================================

def _find_project_root(start: Path) -> Path:
    """
    app.py がある階層、または template/templates ディレクトリがある階層を
    「プロジェクトルート」とみなして探索する。
    """
    cur = start
    for _ in range(6):
        if (cur / "app.py").is_file():
            return cur
        if (cur / "template").is_dir() or (cur / "templates").is_dir():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return start


def _resolve_tool_path(path_hint: str, default_name: str) -> str:
    """
    run.bat / gms_progress.ps1 などのツールパスを解決する。
      1) ユーザー指定 path_hint（存在すれば絶対パスに解決）
      2) プロジェクトルートの template/default_name
      3) プロジェクトルートの templates/default_name
      4) カレントディレクトリの template/templates
    """
    if path_hint:
        p = Path(path_hint).expanduser()
        if p.is_file():
            return str(p.resolve())

    here = Path(__file__).resolve()
    project_root = _find_project_root(here.parent)

    candidates: List[Path] = []
    # project_root 基準
    candidates.append(project_root / "template" / default_name)
    candidates.append(project_root / "templates" / default_name)
    # CWD 基準 (streamlit run app.py のディレクトリ）
    cwd = Path.cwd()
    candidates.append(cwd / "template" / default_name)
    candidates.append(cwd / "templates" / default_name)

    for c in candidates:
        if c.is_file():
            return str(c.resolve())

    msg = f"{default_name} が見つかりませんでした。探索した候補:\n" + "\n".join(
        f"  - {c}" for c in candidates
    )
    raise FileNotFoundError(msg)


# ============================================================
# UMA_summary / post_MM から QM 用構造選抜 & QM フォルダ作成
# ============================================================

def _read_uma_summary(id_dir: Path) -> Optional[pd.DataFrame]:
    """
    ID フォルダ直下の UMA_summary.xlsx を読み込む。
    """
    uma_path = id_dir / "UMA_summary.xlsx"
    if not uma_path.is_file():
        return None
    try:
        return pd.read_excel(uma_path)
    except Exception:
        return None


def _select_qm_candidates_from_uma(df_uma: pd.DataFrame, max_jobs: int) -> List[Dict[str, Any]]:
    """
    UMA_summary の DataFrame から QM に回す候補構造を選ぶ。

    想定カラム:
      - status            : "ok" / "error" / "aborted"
      - xyz_out           : UMA 最適化後 xyz のパス
      - energy_final_eV   : 最終エネルギー（任意）
    """
    cols = set(df_uma.columns)
    if "status" not in cols or "xyz_out" not in cols:
        return []

    df = df_uma.copy()
    df = df[(df["status"] == "ok") & df["xyz_out"].notna()]
    if df.empty:
        return []

    if "energy_final_eV" in cols:
        df = df.sort_values("energy_final_eV", ascending=True)

    cand = df.head(max_jobs)

    candidates: List[Dict[str, Any]] = []
    for _, row in cand.iterrows():
        xyz_path = str(row["xyz_out"])
        if not Path(xyz_path).is_file():
            continue
        candidates.append(
            {
                "xyz_out": xyz_path,
                "energy_final_eV": float(row["energy_final_eV"])
                if "energy_final_eV" in cols and pd.notna(row["energy_final_eV"])
                else None,
            }
        )
    return candidates


def _prepare_qm_from_uma(id_dir: Path, max_jobs: int = 3) -> List[Dict[str, Any]]:
    """
    UMA_summary.xlsx をもとに QM 用の QM/gms1..gmsN フォルダを作成し、
    UMA 最適化後 xyz を geom.xyz としてコピーする。
    """
    df_uma = _read_uma_summary(id_dir)
    if df_uma is None:
        return []

    candidates = _select_qm_candidates_from_uma(df_uma, max_jobs=max_jobs)
    if not candidates:
        return []

    qm_root = id_dir / "QM"
    qm_root.mkdir(parents=True, exist_ok=True)

    selected: List[Dict[str, Any]] = []
    for i, cand in enumerate(candidates, start=1):
        job_dir = qm_root / f"gms{i}"
        job_dir.mkdir(parents=True, exist_ok=True)

        src_xyz = Path(cand["xyz_out"])
        dst_xyz = job_dir / "geom.xyz"

        try:
            shutil.copy2(src_xyz, dst_xyz)
        except Exception:
            continue

        selected.append(
            {
                "job": i,
                "job_dir": str(job_dir),
                "geom_xyz": str(dst_xyz),
                "energy_final_eV": cand.get("energy_final_eV"),
            }
        )

    # ID 直下に QM_summary.xlsx（空）を一旦置いておく（存在確認用）
    if selected:
        try:
            (id_dir / "QM_summary.xlsx").write_bytes(b"")
        except Exception:
            pass

    return selected


def _prepare_qm_from_post_mm(id_dir: Path, max_jobs: int = 3) -> List[Dict[str, Any]]:
    """
    UMA を使わない場合:
      output/<ID>/post_MM/*.xyz を元に QM/gms1..gmsN を作成し、
      conf_0001.xyz などを geom.xyz としてコピーする。
    """
    post_dir = id_dir / "post_MM"
    if not post_dir.is_dir():
        return []

    xyz_files = sorted(post_dir.glob("*.xyz"))
    if not xyz_files:
        return []

    # 最大 max_jobs 個だけ使用（ファイル名順＝エネルギー昇順を想定）
    xyz_files = xyz_files[:max_jobs]

    qm_root = id_dir / "QM"
    qm_root.mkdir(parents=True, exist_ok=True)

    selected: List[Dict[str, Any]] = []
    for i, src in enumerate(xyz_files, start=1):
        job_dir = qm_root / f"gms{i}"
        job_dir.mkdir(parents=True, exist_ok=True)

        dst_xyz = job_dir / "geom.xyz"
        try:
            shutil.copy2(src, dst_xyz)
        except Exception:
            continue

        selected.append(
            {
                "job": i,
                "job_dir": str(job_dir),
                "geom_xyz": str(dst_xyz),
                "energy_final_eV": None,
            }
        )

    if selected:
        try:
            (id_dir / "QM_summary.xlsx").write_bytes(b"")
        except Exception:
            pass

    return selected


# ============================================================
# GAMESS .out 解析用 正規表現など
# ============================================================

_GEOM_HEADER_RE = re.compile(r"^\s*\*{3,}\s*EQUILIBRIUM GEOMETRY LOCATED", re.IGNORECASE)
_COORD_HEADER_RE = re.compile(r"^\s*ATOM\s+CHARGE\s+X\s+Y\s+Z", re.IGNORECASE)
_TOTAL_ENE_RE = re.compile(r"TOTAL ENERGY\s*=\s*([-\d\.]+)", re.IGNORECASE)

# NSTEP 打ち切り（最大ステップ到達）を検知するためのパターン
_MAX_STEP_PATTERNS = [
    r"MAXIMUM NUMBER OF OPTIMIZATION CYCLES REACHED",
    r"MAXIMUM NUMBER OF STEPS TAKEN",
    r"EXCESS NUMBER OF GEOMETRY ITERATIONS",
    r"FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS",
]


def _check_gamess_max_step(out_path: Path) -> bool:
    """
    .out を読んで「最大ステップ数到達」による終了かどうかをざっくり判定する。

    NSTEP=100 に達した場合のみ True を返し、
    他のエラー（SCF 不収束など）の場合は False のままにして、従来通りエラーとする。
    """
    if not out_path.is_file():
        return False
    text = out_path.read_text(errors="ignore")
    for pat in _MAX_STEP_PATTERNS:
        if re.search(pat, text, re.IGNORECASE):
            return True
    return False


def _parse_equilibrium_geometry_and_energy(
    out_path: Path,
) -> Tuple[Optional[float], List[Tuple[str, float, float, float]]]:
    """
    GAMESS の .out から
      - 'EQUILIBRIUM GEOMETRY LOCATED' ブロックの座標
      - TOTAL ENERGY
    を取り出す。
    """
    if not out_path.is_file():
        return None, []

    text = out_path.read_text(errors="ignore")
    lines = text.splitlines()

    # --- エネルギー（ファイル全体から最後の TOTAL ENERGY を拾う） ---
    energies = [float(m.group(1)) for m in _TOTAL_ENE_RE.finditer(text)]
    energy_h: Optional[float] = energies[-1] if energies else None

    # --- ジオメトリ: 最後の "EQUILIBRIUM GEOMETRY LOCATED" を探す ---
    geom_start: Optional[int] = None
    for i in range(len(lines) - 1, -1, -1):
        if _GEOM_HEADER_RE.search(lines[i]):
            geom_start = i
            break
    if geom_start is None:
        return energy_h, []

    # "ATOM  CHARGE  X  Y  Z" ヘッダを探す（少し広めに 200 行）
    coord_header: Optional[int] = None
    for j in range(geom_start, min(geom_start + 200, len(lines))):
        if _COORD_HEADER_RE.search(lines[j]):
            coord_header = j
            break
    if coord_header is None:
        return energy_h, []

    # ヘッダの次行以降で、最初の "-----" ラインを探してその次を座標開始とする
    idx = coord_header + 1
    while idx < len(lines) and "-" not in lines[idx]:
        idx += 1
    if idx >= len(lines):
        return energy_h, []
    idx += 1  # ダッシュ行の次が最初の座標行

    coords: List[Tuple[str, float, float, float]] = []

    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            break

        parts = line.split()
        if len(parts) < 4:
            break

        sym = parts[0]

        # 行末 3 つを X, Y, Z とみなす
        try:
            x = float(parts[-3])
            y = float(parts[-2])
            z = float(parts[-1])
        except Exception:
            break

        coords.append((sym, x, y, z))
        idx += 1

    return energy_h, coords


def _parse_last_geometry(out_path: Path) -> List[Tuple[str, float, float, float]]:
    """
    最後の最適化ステップの座標を抜き出す。
    EQUILIBRIUM GEOMETRY LOCATED が無い（=未収束）場合のリスタート用。

    GAMESS .out 内の
      'COORDINATES OF ALL ATOMS ARE (ANGS)'
    ブロックを後ろから探し、そのテーブルを座標として読み取る。
    """
    if not out_path.is_file():
        return []

    text = out_path.read_text(errors="ignore")
    lines = text.splitlines()

    header_re = re.compile(r"COORDINATES OF ALL ATOMS ARE\s*\(ANGS\)", re.IGNORECASE)

    header_idx: Optional[int] = None
    for i in range(len(lines) - 1, -1, -1):
        if header_re.search(lines[i]):
            header_idx = i
            break

    if header_idx is None:
        return []

    idx = header_idx + 3
    coords: List[Tuple[str, float, float, float]] = []

    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            break

        parts = line.split()
        # 例: "1  C   6.0   x   y   z" → ["1","C","6.0","x","y","z"]
        if len(parts) < 6:
            break

        sym = parts[1]
        try:
            x = float(parts[-3])
            y = float(parts[-2])
            z = float(parts[-1])
        except Exception:
            break

        coords.append((sym, x, y, z))
        idx += 1

    return coords


def _write_xyz_from_coords(
    coords: List[Tuple[str, float, float, float]],
    energy_h: Optional[float],
    out_xyz: Path,
    comment: str,
) -> None:
    """
    座標とエネルギーからシンプルな XYZ ファイルを書き出す。
    """
    n = len(coords)
    out_lines: List[str] = []
    out_lines.append(str(n))
    if energy_h is not None:
        out_lines.append(f"{comment}  E={energy_h:.10f} Hartree")
    else:
        out_lines.append(comment)
    for sym, x, y, z in coords:
        out_lines.append(f"{sym:2s}  {x: .10f}  {y: .10f}  {z: .10f}")
    out_xyz.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


# ============================================================
# サブプロセス実行 + 即時中断監視
# ============================================================

def _run_subprocess_with_abort(
    cmd: List[str],
    cwd: Path,
    runner_log: Path,
    slot: str,
    attempt: int,
    is_abort_now: Callable[[], bool],
) -> int:
    """
    GAMESS 実行用のサブプロセスを起動し、abort フラグを 1 秒ごとに監視する。
    中断が指示された場合はプロセスを強制終了して RuntimeError を投げる。
    正常終了した場合は returncode を返す。
    """
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    stdout_accum = ""
    stderr_accum = ""

    try:
        while True:
            try:
                out, err = proc.communicate(timeout=1.0)
                stdout_accum += out or ""
                stderr_accum += err or ""
                break
            except subprocess.TimeoutExpired:
                if is_abort_now():
                    # 通常の terminate
                    try:
                        proc.terminate()
                    except Exception:
                        pass
                    # Windows では子プロセスごと kill（失敗しても無視）
                    try:
                        subprocess.run(
                            ["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                        )
                    except Exception:
                        pass
                    # ログに中断を記録
                    try:
                        with runner_log.open("a", encoding="utf-8") as f:
                            f.write(
                                f"[{slot}] attempt={attempt} aborted by user via abort flag.\n"
                            )
                    except Exception:
                        pass
                    raise RuntimeError(
                        f"GAMESS job '{slot}' aborted by user (attempt={attempt})."
                    )
    finally:
        returncode = proc.returncode

        # ログ出力
        try:
            with runner_log.open("a", encoding="utf-8") as f:
                f.write(f"[{slot}] attempt={attempt} CMD: {' '.join(cmd)}\n")
                f.write(f"[{slot}] returncode: {returncode}\n")
                if stdout_accum:
                    f.write(f"[{slot}] STDOUT:\n{stdout_accum}\n")
                if stderr_accum:
                    f.write(f"[{slot}] STDERR:\n{stderr_accum}\n")
                f.write("\n")
        except Exception:
            pass

    return returncode


# ============================================================
# GAMESS 実行（run.bat 呼び出し） + NSTEP リスタート + 即時中断
# ============================================================

def _run_gamess_jobs_sequential(
    id_dir: Path,
    qm_cfg: Any,
    slots: Sequence[str],
    is_abort_now: Optional[Callable[[], bool]] = None,
) -> None:
    """
    QM/gmsX/slot.inp を run.bat で 1 つずつ直列実行する。

    - run.bat / gms_progress.ps1 のパスを解決
    - 実行コマンドと returncode を logs/qm_runner.log に記録
    - 各スロットについて:
        1. 一度実行
        2. EQUILIBRIUM GEOMETRY LOCATED が見つかれば収束とみなして次のスロットへ
        3. 見つからず、かつ .out に NSTEP による打切りメッセージがある場合のみ
           最後の座標ブロックから geom.xyz を作り直し、入力を再生成して再実行
        4. NSTEP 以外の理由で未収束なら即エラー
    さらに、実行中は 1 秒ごとに abort フラグを監視し、中断が指示されたら
    現在実行中の GAMESS を kill して RuntimeError を投げる。
    """
    if is_abort_now is None:
        def is_abort_now() -> bool:
            return False

    qm_root = id_dir / "QM"
    logs_dir = id_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    runner_log = logs_dir / "qm_runner.log"

    # run.bat / gms_progress.ps1 のパス解決（絶対パス）
    bat = _resolve_tool_path(getattr(qm_cfg, "run_bat_path", "") or "", "run.bat")
    ps1 = _resolve_tool_path(getattr(qm_cfg, "ps1_path", "") or "", "gms_progress.ps1")

    ncpus = int(getattr(qm_cfg, "ncpus", 4))

    # NSTEP 打ち切りからの自動リスタート最大回数
    max_restarts = 2  # 初回 + リスタート2回 = 最大3実行

    # gms.inp 再生成用パラメータ（リスタート時にも使う）
    tpl_hint = getattr(qm_cfg, "template_path", "") or ""
    smiles_comment = getattr(qm_cfg, "smiles_comment", "")
    pointgroup = getattr(qm_cfg, "pointgroup", "C1")
    icharge = int(getattr(qm_cfg, "icharge", 0))
    mult = int(getattr(qm_cfg, "mult", 1))

    for slot in slots:
        # スロット開始前にも abort をチェック
        if is_abort_now():
            raise RuntimeError(
                f"GAMESS jobs aborted by user before starting slot '{slot}'."
            )

        sdir = qm_root / slot
        out_file = sdir / f"{slot}.out"

        for attempt in range(max_restarts + 1):
            inp = sdir / f"{slot}.inp"
            if not inp.is_file():
                # .inp が無ければこのスロットは諦める
                raise RuntimeError(
                    f"GAMESS job '{slot}' の入力ファイルが見つかりません: {inp}"
                )

            cmd = [
                "cmd",
                "/c",
                bat,
                str(inp.resolve()),
                str(ncpus),
                str(logs_dir.resolve()),
                str(ps1),
            ]

            # abort 監視付きで実行
            returncode = _run_subprocess_with_abort(
                cmd=cmd,
                cwd=sdir,
                runner_log=runner_log,
                slot=slot,
                attempt=attempt,
                is_abort_now=is_abort_now,
            )

            # RETURN CODE/OUT ファイルの存在チェック
            if returncode != 0:
                raise RuntimeError(
                    f"GAMESS job '{slot}' failed (attempt={attempt}, returncode={returncode}).\n"
                    f"see {runner_log} for details."
                )

            if not out_file.is_file():
                raise RuntimeError(
                    f"GAMESS job '{slot}' finished with returncode 0, "
                    f"but output file not found: {out_file}"
                )

            # --- ここから収束判定と NSTEP リスタート処理 ---
            energy_h, coords_eq = _parse_equilibrium_geometry_and_energy(out_file)
            if coords_eq:
                # EQUILIBRIUM GEOMETRY LOCATED を取得できた → 収束
                break

            # ここに来るのは「EQUILIBRIUM GEOMETRY LOCATED が無い」ケース。
            # NSTEP 打ち切りかどうかを判定する。
            if not _check_gamess_max_step(out_file):
                # NSTEP 以外の理由で未収束 → 即エラー
                raise RuntimeError(
                    f"GAMESS job '{slot}' did not converge and it is not a max-step (NSTEP) stop.\n"
                    f"see {out_file} for details."
                )

            # ここに来るのは「NSTEP に達して止まった」ケース。
            if attempt >= max_restarts:
                # 既に最大リスタート回数を超えている → エラーであきらめる
                raise RuntimeError(
                    f"GAMESS job '{slot}' reached max NSTEP and did not converge "
                    f"after {max_restarts + 1} attempts."
                )

            # リスタート: 最後の座標ブロックを取り出して geom.xyz を作り直し、
            # それを用いた新しい .inp を生成する。
            last_coords = _parse_last_geometry(out_file)
            if not last_coords:
                raise RuntimeError(
                    f"GAMESS job '{slot}' reached max NSTEP but last geometry "
                    f"could not be parsed from {out_file}."
                )

            # geom.xyz を上書き
            geom_xyz = sdir / "geom.xyz"
            _write_xyz_from_coords(
                last_coords,
                energy_h=None,
                out_xyz=geom_xyz,
                comment=f"restart {slot} attempt={attempt+1}",
            )

            # 対象スロットだけ gms.inp を再生成
            make_inp_for_slots(
                id_dir=str(id_dir),
                gms_template_path=tpl_hint,
                smiles_comment=smiles_comment or slot,
                pointgroup=pointgroup,
                icharge=icharge,
                mult=mult,
                slots=[slot],
            )

        # for attempt が正常に break した場合のみ次の slot へ


# ============================================================
# QM 結果のまとめ & result/ へのコピー
# ============================================================

def _postprocess_qm_for_id(
    mol_id: str,
    id_dir: Path,
    slots: Sequence[str],
    result_root: Path,
) -> QMResult:
    """
    QM/gmsX/*.out から最終構造とエネルギーを抜き出し、
      - QM_summary.xlsx（ID 直下）
      - gmsX_final.xyz（各スロット）
      - result/<ID>/<ID>_best.xyz（最安定構造）
    を生成する。
    """
    qm_root = id_dir / "QM"
    summary_records: List[Dict[str, Any]] = []

    for slot in slots:
        sdir = qm_root / slot
        out_file = sdir / f"{slot}.out"
        energy_h, coords = _parse_equilibrium_geometry_and_energy(out_file)

        final_xyz_path: Optional[str] = None
        if coords:
            final_xyz = sdir / f"{slot}_final.xyz"
            _write_xyz_from_coords(coords, energy_h, final_xyz, comment=f"{mol_id} {slot}")
            final_xyz_path = str(final_xyz)

        summary_records.append(
            {
                "slot": slot,
                "out_path": str(out_file),
                "energy_hartree": energy_h,
                "final_xyz": final_xyz_path or "",
            }
        )

    # QM_summary.xlsx
    summary_path = id_dir / "QM_summary.xlsx"
    if summary_records:
        df = pd.DataFrame(summary_records)
    else:
        df = pd.DataFrame(columns=["slot", "out_path", "energy_hartree", "final_xyz"])
    df.to_excel(summary_path, index=False, sheet_name="QM")

    # 最安定構造を result/ にコピー
    best_energy: Optional[float] = None
    best_xyz: Optional[str] = None
    for rec in summary_records:
        e = rec["energy_hartree"]
        xyz = rec["final_xyz"]
        if e is None or not xyz:
            continue
        if best_energy is None or e < best_energy:
            best_energy = e
            best_xyz = xyz

    best_xyz_out: Optional[str] = None
    if best_xyz:
        src = Path(best_xyz)
        result_id_dir = result_root / mol_id
        result_id_dir.mkdir(parents=True, exist_ok=True)
        dst = result_id_dir / f"{mol_id}_best.xyz"
        shutil.copy2(src, dst)
        best_xyz_out = str(dst)

    return QMResult(
        id=mol_id,
        n_jobs=len(slots),
        summary_path=str(summary_path),
        best_energy_h=best_energy,
        best_xyz_path=best_xyz_out,
    )


# ============================================================
# 公開関数: run_qm_for_job
# ============================================================

def run_qm_for_job(
    mol_id: str,
    uma_summary: Union[str, UMAResult, None],   # UMA あり: str/UMAResult, UMA なし: None
    qm_cfg: Any,
    out_cfg: Any,
    is_abort_now: Optional[Callable[[], bool]] = None,
) -> QMResult:
    """
    app.py から呼ばれるエントリポイント。

    uma_summary:
        - UMA_summary.xlsx のパス（str）
        - UMAResult オブジェクト
        - UMA を使わない場合は None

    フロー:
      UMA あり:
        1. UMA_summary.xlsx から QM/gms1..gmsN を構築
      UMA なし:
        1. output/<ID>/post_MM/*.xyz から QM/gms1..gmsN を構築

      2. make_inp_for_slots() で gms.inp テンプレートから .inp を生成
      3. run.bat を用いて gmsX を直列実行（NSTEP 打ち切り時は自動リスタート + 即時中断監視）
      4. .out から最終構造/エネルギーを抽出
      5. QM_summary.xlsx と result/<ID>/<ID>_best.xyz を生成
    """
    base_dir = getattr(out_cfg, "base_dir", "./output")

    # ---- UMA あり or なしで分岐 ----
    if uma_summary is None:
        # UMA を使わないケース: output/<ID>/post_MM から候補を取る
        id_dir = (Path(base_dir) / mol_id).resolve()
        selected = _prepare_qm_from_post_mm(id_dir, max_jobs=3)
        source_kind = "post_MM"
    else:
        # UMA あり
        if isinstance(uma_summary, UMAResult):
            summary_path_str = getattr(uma_summary, "summary_path", "") or ""
            if summary_path_str:
                uma_path = Path(summary_path_str)
            else:
                uma_path = Path(base_dir) / mol_id / "UMA_summary.xlsx"
        else:
            uma_path = Path(str(uma_summary))

        if uma_path.is_file():
            id_dir = uma_path.parent.resolve()
        else:
            id_dir = (Path(base_dir) / mol_id).resolve()

        selected = _prepare_qm_from_uma(id_dir, max_jobs=3)
        source_kind = "UMA"

    if not selected:
        if source_kind == "UMA":
            raise RuntimeError(f"[{mol_id}] UMA 正常終了構造が無いため QM 構造準備をスキップします。")
        else:
            raise RuntimeError(
                f"[{mol_id}] post_MM から QM 用構造を準備できませんでした。"
                f" output/{mol_id}/post_MM を確認してください。"
            )

    slots = [f"gms{rec['job']}" for rec in selected]

    # gms.inp 生成
    tpl_hint = getattr(qm_cfg, "template_path", "") or ""
    smiles_comment = getattr(qm_cfg, "smiles_comment", mol_id)
    pointgroup = getattr(qm_cfg, "pointgroup", "C1")
    icharge = int(getattr(qm_cfg, "icharge", 0))
    mult = int(getattr(qm_cfg, "mult", 1))

    make_inp_for_slots(
        id_dir=str(id_dir),
        gms_template_path=tpl_hint,
        smiles_comment=smiles_comment,
        pointgroup=pointgroup,
        icharge=icharge,
        mult=mult,
        slots=slots,
    )

    # GAMESS を直列実行（NSTEP 打ち切り時は自動リスタート + 即時中断監視）
    _run_gamess_jobs_sequential(id_dir, qm_cfg, slots, is_abort_now=is_abort_now)

    # post 処理: QM_summary.xlsx & result/<ID> 出力
    result_root = Path.cwd() / "result"
    qm_result = _postprocess_qm_for_id(mol_id, id_dir, slots, result_root)

    return qm_result
