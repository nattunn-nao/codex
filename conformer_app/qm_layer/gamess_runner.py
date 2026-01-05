from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Sequence

from conformer_app.core.config import QMConfig


# ======================================
# パス解決ユーティリティ
# ======================================


def _project_root() -> Path:
    """
    re_production/ をプロジェクトルートとみなす。
    （conformer_app/qm_layer/ から 2 つ上の階層）
    """
    return Path(__file__).resolve().parents[2]


def _resolve_tool_path(path_hint: str | None, default_name: str) -> Path:
    """
    run.bat / gms_progress.ps1 の場所を解決する。

    優先順位:
      1) QMConfig から渡された path_hint
         - 相対パスならカレントディレクトリ基準 (re_production/) で解決
      2) <プロジェクトルート>/template/<default_name>
      3) CWD/template/<default_name>
    """
    candidates: list[Path] = []

    if path_hint:
        p = Path(path_hint)
        if not p.is_absolute():
            p = Path.cwd() / p
        candidates.append(p)

    root = _project_root()
    candidates.append(root / "template" / default_name)
    candidates.append(Path.cwd() / "template" / default_name)

    for p in candidates:
        if p.exists():
            return p

    raise FileNotFoundError(
        f"{default_name} が見つかりませんでした: "
        f"hint={path_hint!r}, tried={[str(c) for c in candidates]}"
    )


# ======================================
# 単一ジョブ実行
# ======================================


def _run_single_gamess(
    inp: Path,
    run_bat: Path,
    ps1: Path,
    log_dir: Path,
    ncpus: int,
) -> dict:
    """
    1 つの .inp に対して GAMESS を実行し、結果を辞書で返す。
    """
    inp = inp.resolve()
    job_dir = inp.parent
    job_name = inp.stem
    log_dir = log_dir.resolve()

    cmd = [
        "cmd",
        "/c",
        str(run_bat),
        str(inp),
        str(ncpus),
        str(log_dir),
        str(ps1.resolve()),
    ]

    proc = subprocess.run(
        cmd,
        cwd=str(job_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    out_path = inp.with_suffix(".out")

    return {
        "name": job_name,
        "cmd": cmd,
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "out_path": out_path,
    }


# ======================================
# メイン: GAMESS ジョブ並列実行
# ======================================


def run_gamess_jobs(
    id_dir: Path,
    inp_paths: Sequence[Path],
    qm_cfg: QMConfig,
) -> List[Path]:
    """
    1つの ID について、用意された複数の .inp を
    **最大 max_parallel 並列**で GAMESS 実行する。

    Parameters
    ----------
    id_dir : Path
        その ID のルートフォルダ（例: ./output/111-200）
    inp_paths : Sequence[Path]
        QM/gms1/gms1.inp, QM/gms2/gms2.inp ... のフルパス
    qm_cfg : QMConfig
        run.bat / ps1 のパス、NCPUS、max_parallel 等を含む設定

    Returns
    -------
    out_paths : list[Path]
        正常終了したジョブの .out パス一覧。
        1つでも失敗ジョブがあれば RuntimeError を投げる。
    """
    id_dir = Path(id_dir)
    inp_list = [Path(p) for p in inp_paths if Path(p).exists()]
    if not inp_list:
        return []

    # ログディレクトリ & ランチャーログ
    log_dir = id_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    runner_log = log_dir / "qm_runner.log"

    # ツールパス解決
    run_bat = _resolve_tool_path(getattr(qm_cfg, "run_bat_path", ""), "run.bat")
    ps1 = _resolve_tool_path(getattr(qm_cfg, "ps1_path", ""), "gms_progress.ps1")

    ncpus = int(getattr(qm_cfg, "ncpus", 4))
    max_parallel = int(getattr(qm_cfg, "max_parallel", 3))
    if max_parallel < 1:
        max_parallel = 1

    results: list[dict] = []

    # 並列実行（ThreadPoolExecutor でシンプルに）
    with ThreadPoolExecutor(max_workers=max_parallel) as ex:
        future_to_inp = {
            ex.submit(_run_single_gamess, inp, run_bat, ps1, log_dir, ncpus): inp
            for inp in inp_list
        }
        for fut in as_completed(future_to_inp):
            res = fut.result()
            results.append(res)

    # ログ書き込み & 成功判定
    out_paths: List[Path] = []
    failed: List[dict] = []

    with runner_log.open("a", encoding="utf-8", errors="ignore") as lf:
        for res in results:
            name = res["name"]
            cmd_str = " ".join(res["cmd"])
            rc = res["returncode"]
            stdout = res.get("stdout") or ""
            stderr = res.get("stderr") or ""
            out_path: Path = res["out_path"]

            lf.write(f"[{name}] CMD: {cmd_str}\n")
            lf.write(f"[{name}] returncode: {rc}\n")
            if stdout.strip():
                lf.write(f"[{name}] STDOUT:\n{stdout}\n")
            if stderr.strip():
                lf.write(f"[{name}] STDERR:\n{stderr}\n")

            if rc == 0 and out_path.exists():
                out_paths.append(out_path)
            else:
                failed.append(res)

    if failed:
        first = failed[0]
        raise RuntimeError(
            f"GAMESS job '{first['name']}' failed (returncode={first['returncode']}). "
            f"see {runner_log} for details."
        )

    return out_paths


def launch_gamess_parallel(
    id_dir: str | Path,
    inp_paths: Sequence[str | Path],
    qm_cfg: QMConfig,
) -> List[Path]:
    """
    旧 confapp.qm_gamess.launch_gamess_parallel 互換のラッパー。
    """
    return run_gamess_jobs(Path(id_dir), [Path(p) for p in inp_paths], qm_cfg)
