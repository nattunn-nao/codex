from __future__ import annotations

import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS
from conformer_app.core.config import UMAConfig, OutputConfig
from conformer_app.uma_layer.calculator import create_uma_calculator
from conformer_app.uma_layer.summary import write_uma_summary

# Windows + fairchem + OpenMP 対策（多重ロードを許容）
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")


# ==========================
# 公開: UMAResult
# ==========================


@dataclass
class UMAResult:
    """
    UMA ステージ全体の結果まとめ。

    Attributes
    ----------
    status : str
        "ok" / "error" / "aborted"
    message : str
        補足メッセージ（エラーや警告など）
    uma_summary_path : str
        <ID>/UMA_summary.xlsx のパス（失敗時は空文字）
    uma_dir : str
        <ID>/UMA/ ディレクトリのパス
    n_structures : int
        UMA 最適化を試みた構造数
    """

    status: str
    message: str = ""
    uma_summary_path: str = ""
    uma_dir: str = ""
    n_structures: int = 0


# ==========================
# 内部ユーティリティ
# ==========================


def _log_side(msg: str) -> None:
    """Streamlit があればサイドバーへ軽くメッセージを出す。"""
    try:
        import streamlit as st  # type: ignore

        st.sidebar.write(msg)
    except Exception:
        pass


def _log_uma_line(log_path: Optional[str], msg: str) -> None:
    """
    UMA 用のログ出力（1 conf ごとに専用の log ファイル）。
    タイムスタンプは付けず、msg をそのまま 1 行書く。
    """
    if not log_path:
        return
    try:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(msg + "\n")
    except Exception:
        # ログ出力失敗は致命的ではないので無視
        pass


def _rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """2つの座標配列の RMSD (Å) を計算する。shape = (N, 3) を想定。"""
    diff = a - b
    return float(np.sqrt((diff * diff).sum() / diff.shape[0]))


def _force_stats(forces: np.ndarray) -> Tuple[float, float]:
    """
    力ベクトル forces (N,3) から f_max, f_rms を計算する。
    f_max : 各原子の |F| の最大値
    f_rms : 各原子の |F| の RMS
    """
    if forces.size == 0:
        return float("inf"), float("inf")
    norms = np.linalg.norm(forces, axis=1)
    f_max = float(norms.max())
    f_rms = float(np.sqrt((norms * norms).mean()))
    return f_max, f_rms


# ==========================
# UMA バッチ最適化本体
# ==========================


def _uma_optimize_xyz_batch_fairchem(
    xyz_files: List[str],
    out_dir: str,
    cfg: UMAConfig,
    id_dir: Optional[str] = None,
    is_abort_now: Optional[Callable[[], bool]] = None,
) -> List[Dict[str, Any]]:
    """
    UMA ポテンシャル + LBFGS で xyz 群を最適化するバッチ処理。
    各構造ごとに専用 log を出力し、step ごとの Energy / f_max / f_rms と
    start_time / end_time / wall_time を記録する。
    """
    if is_abort_now is None:

        def is_abort_now() -> bool:
            return False

    os.makedirs(out_dir, exist_ok=True)

    # cfg.model_name / cfg.model_tag のどちらでも動くようにする
    model_name = getattr(
        cfg, "model_name", getattr(cfg, "model_tag", "uma-s-1p1")
    )
    device = getattr(cfg, "device", "cpu")
    task_name = getattr(cfg, "task_name", "omol")
    fmax_thresh = float(getattr(cfg, "fmax", 0.01))
    max_steps = int(getattr(cfg, "max_steps", 500))
    debug_flag = bool(getattr(cfg, "debug", False))

    # ログ出力先ディレクトリ（IDごとの logs/）
    if id_dir:
        logs_dir = os.path.join(id_dir, "logs")
    else:
        logs_dir = out_dir
    os.makedirs(logs_dir, exist_ok=True)

    batch_start = time.time()

    # Calculator 構築
    calc = create_uma_calculator(
        UMAConfig(
            model_name=model_name,
            device=device,
            task_name=task_name,
        )
    )

    results: List[Dict[str, Any]] = []

    for idx, xyz_path in enumerate(xyz_files, start=1):
        if is_abort_now():
            # この xyz の計算すら始めずに中断
            base = os.path.basename(xyz_path)
            base_noext, _ = os.path.splitext(base)
            log_path = os.path.join(logs_dir, f"UMA_{base_noext}.log")
            _log_uma_line(log_path, f"# UMA aborted before start idx={idx} file={base}")
            results.append(
                {
                    "idx": idx,
                    "xyz_in": xyz_path,
                    "xyz_out": "",
                    "status": "aborted",
                    "energy_init_eV": None,
                    "energy_final_eV": None,
                    "rmsd_final_A": None,
                    "message": "aborted before start",
                }
            )
            break

        base = os.path.basename(xyz_path)
        base_noext, _ = os.path.splitext(base)
        log_path = os.path.join(logs_dir, f"UMA_{base_noext}.log")

        file_start = time.time()
        _log_uma_line(log_path, "")
        _log_uma_line(log_path, f"# ==== UMA start idx={idx} file={base} ====")
        _log_uma_line(
            log_path,
            f"# model={model_name}, device={device}, fmax={fmax_thresh}, "
            f"max_steps={max_steps}, debug={debug_flag}",
        )
        _log_uma_line(
            log_path,
            f"# start_time {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(file_start))}",
        )
        _log_uma_line(log_path, "# step  energy(eV)           f_max(eV/A)       f_rms(eV/A)")

        try:
            atoms = read(xyz_path)
            # charge / spin : UMAConfig から優先的に設定し、無ければデフォルト
            atoms.info = atoms.info or {}
            atoms.info.setdefault("charge", getattr(cfg, "charge", 0))
            atoms.info.setdefault("spin", getattr(cfg, "spin", 1))

            atoms.calc = calc

            # 初期座標
            pos0 = atoms.get_positions().copy()

            opt = LBFGS(atoms, logfile=None)
            step = 0
            converged = False
            e0 = None
            e_last = None

            while True:
                if is_abort_now():
                    raise KeyboardInterrupt("UMA optimization aborted by user")

                # 現在の状態で energy / forces を評価
                e_curr = float(atoms.get_potential_energy())
                forces = atoms.get_forces()
                f_max, f_rms = _force_stats(forces)

                step += 1
                if e0 is None:
                    e0 = e_curr
                e_last = e_curr

                # ログ出力：step, energy, f_max, f_rms（すべて実数表記）
                _log_uma_line(
                    log_path,
                    f"{step:5d}  {e_curr: .10f}  {f_max: .8f}  {f_rms: .8f}",
                )

                # 収束判定
                if f_max < fmax_thresh or step >= max_steps:
                    converged = f_max < fmax_thresh
                    break

                # 1 step 進める
                opt.step()

            # 最終エネルギーと座標
            e1 = e_last if e_last is not None else float(atoms.get_potential_energy())
            pos1 = atoms.get_positions().copy()

            # RMSD (Å)
            try:
                rmsd_val = _rmsd(pos0, pos1)
            except Exception:
                rmsd_val = None

            # 出力 xyz
            out_xyz = os.path.join(out_dir, base)
            write(out_xyz, atoms)

            file_end = time.time()
            wall = file_end - file_start

            _log_uma_line(
                log_path,
                f"# end_time   {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(file_end))}",
            )
            _log_uma_line(log_path, f"# wall_time  {wall:.3f} s")
            _log_uma_line(
                log_path,
                f"# converged={converged}, steps={step}, "
                f"E0={e0 if e0 is not None else 'NA'} eV, E1={e1:.10f} eV, "
                f"RMSD={rmsd_val if rmsd_val is not None else 'NA'} Å",
            )
            _log_uma_line(log_path, f"# ==== UMA end idx={idx} file={base} ====")

            results.append(
                {
                    "idx": idx,
                    "xyz_in": xyz_path,
                    "xyz_out": out_xyz,
                    "status": "ok",
                    "energy_init_eV": float(e0) if e0 is not None else None,
                    "energy_final_eV": float(e1),
                    "rmsd_final_A": rmsd_val,
                    "message": "",
                }
            )

        except KeyboardInterrupt:
            file_end = time.time()
            wall = file_end - file_start
            _log_uma_line(log_path, f"# UMA aborted during optimization idx={idx} file={base}")
            _log_uma_line(
                log_path,
                f"# end_time   {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(file_end))}",
            )
            _log_uma_line(log_path, f"# wall_time  {wall:.3f} s")
            results.append(
                {
                    "idx": idx,
                    "xyz_in": xyz_path,
                    "xyz_out": "",
                    "status": "aborted",
                    "energy_init_eV": None,
                    "energy_final_eV": None,
                    "rmsd_final_A": None,
                    "message": "aborted during optimization",
                }
            )
            break

        except Exception as e:
            file_end = time.time()
            wall = file_end - file_start
            _log_uma_line(log_path, f"# UMA ERROR idx={idx} file={base}: {repr(e)}")
            _log_uma_line(
                log_path,
                f"# end_time   {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(file_end))}",
            )
            _log_uma_line(log_path, f"# wall_time  {wall:.3f} s")
            results.append(
                {
                    "idx": idx,
                    "xyz_in": xyz_path,
                    "xyz_out": "",
                    "status": "error",
                    "energy_init_eV": None,
                    "energy_final_eV": None,
                    "rmsd_final_A": None,
                    "message": repr(e),
                }
            )

    batch_end = time.time()
    batch_wall = batch_end - batch_start
    # バッチ全体の情報は最後に summary として 1 行だけ別ファイルに出しておく
    if id_dir:
        logs_dir = os.path.join(id_dir, "logs")
        batch_log = os.path.join(logs_dir, "UMA_batch.log")
        _log_uma_line(
            batch_log,
            f"UMA batch end: n_results={len(results)}, wall_time={batch_wall:.3f} s",
        )

    return results


# ==========================
# 公開関数: run_uma_for_job
# ==========================


def run_uma_for_job(
    mol_id: str,
    mm_result: Any,
    uma_cfg: UMAConfig,
    out_cfg: Optional[OutputConfig] = None,
    is_abort_now: Optional[Callable[[], bool]] = None,
) -> UMAResult:
    """
    1 つの ID について、MM 出力 xyz 群を UMA ポテンシャルで最適化し、
    UMA_summary.xlsx を作成する高レベル関数。

    互換性ポリシー:
      - まず mm_result.xyz_files があればそれを優先して UMA に渡す
        （= post_MM/*.xyz を想定）
      - 無ければ output/<ID>/post_MM/*.xyz を探す
      - それでも無ければ従来どおり output/<ID>/MM/*.xyz を見る
    """
    # ---- base_dir / id_dir の決定（out_cfg が無くても動くように）----
    if out_cfg is None:
        base_dir = Path("./output")
    else:
        # OutputConfig.base_dir を優先し、無ければ文字列として解釈
        base_dir_val = getattr(out_cfg, "base_dir", None)
        if base_dir_val:
            base_dir = Path(base_dir_val)
        else:
            base_dir = Path(str(out_cfg))

    id_dir = base_dir / mol_id
    mm_dir = id_dir / "MM"
    post_mm_dir = id_dir / "post_MM"
    uma_dir = id_dir / "UMA"

    # === UMA に渡す xyz ファイルの決定ロジック ===
    xyz_files: List[str] = []

    # 1) 新しい MMResult 形式：mm_result.xyz_files（post_MM/*.xyz を想定）
    if hasattr(mm_result, "xyz_files") and mm_result.xyz_files:
        xyz_files = list(mm_result.xyz_files)

    # 2) mm_result に情報が無い場合 → post_MM を優先して使う
    if not xyz_files and post_mm_dir.is_dir():
        xyz_files = sorted(str(p) for p in post_mm_dir.glob("*.xyz"))

    # 3) それでも空なら最後の fallback として MM/*.xyz（旧挙動）
    if not xyz_files and mm_dir.is_dir():
        xyz_files = sorted(str(p) for p in mm_dir.glob("*.xyz"))

    if not xyz_files:
        msg = (
            f"[{mol_id}] UMA に渡す xyz が見つからないため UMA をスキップします "
            f"(post_MM={post_mm_dir}, MM={mm_dir})."
        )
        _log_side(msg)
        return UMAResult(
            status="error",
            message=msg,
            uma_summary_path="",
            uma_dir=str(uma_dir),
            n_structures=0,
        )

    try:
        results = _uma_optimize_xyz_batch_fairchem(
            xyz_files=xyz_files,
            out_dir=str(uma_dir),
            cfg=uma_cfg,
            id_dir=str(id_dir),
            is_abort_now=is_abort_now,
        )
    except KeyboardInterrupt:
        msg = f"[{mol_id}] UMA がユーザー操作により中断されました。"
        _log_side(msg)
        return UMAResult(
            status="aborted",
            message=msg,
            uma_summary_path="",
            uma_dir=str(uma_dir),
            n_structures=0,
        )
    except Exception as e:
        msg = f"[{mol_id}] UMA 実行中に例外発生: {e}"
        _log_side(msg)
        return UMAResult(
            status="error",
            message=msg,
            uma_summary_path="",
            uma_dir=str(uma_dir),
            n_structures=0,
        )

    summary_path = write_uma_summary(str(id_dir), results)

    # 結果のステータス集計
    any_ok = any(r.get("status") == "ok" for r in results)
    any_abort = any(r.get("status") == "aborted" for r in results)

    if any_ok and not any_abort:
        status = "ok"
        message = ""
    elif any_ok and any_abort:
        status = "ok"
        message = "一部の構造が UMA 中断/失敗しましたが、正常終了した構造もあります。"
    elif any_abort and not any_ok:
        status = "aborted"
        message = "すべての UMA 計算が中断されました。"
    else:
        status = "error"
        message = "UMA 計算がすべて失敗しました。"

    return UMAResult(
        status=status,
        message=message,
        uma_summary_path=summary_path or "",
        uma_dir=str(uma_dir),
        n_structures=len(results),
    )
