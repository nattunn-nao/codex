from __future__ import annotations

import os
import time
from typing import Any, Callable, Dict, List, Optional

import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator


# ==========================
# 内部ユーティリティ
# ==========================


def _build_calculator(model_name: str, device: str, task_name: str) -> FAIRChemCalculator:
    """
    fairchem の事前学習モデルから Calculator を構築する。
    """
    predictor = pretrained_mlip.get_predict_unit(model_name, device=device)
    calc = FAIRChemCalculator(predictor, task_name=task_name)
    return calc


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


def _force_stats(forces: np.ndarray) -> tuple[float, float]:
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
# 公開関数: UMA バッチ最適化
# ==========================


def uma_optimize_xyz_batch_fairchem(
    xyz_files: List[str],
    out_dir: str,
    cfg,
    id_dir: Optional[str] = None,
    is_abort_now: Optional[Callable[[], bool]] = None,
    chunk_steps: int = 20,  # 互換性のため残すが、実装上は各 step ごとにチェック
) -> List[Dict[str, Any]]:
    """
    UMA ポテンシャル + LBFGS で xyz 群を最適化するバッチ処理。
    各構造ごとに専用 log を出力し、step ごとの Energy / f_max / f_rms と
    start_time / end_time / wall_time を記録する。

    Parameters
    ----------
    xyz_files : list of str
        入力 xyz ファイル群（必ず呼び出し側で post_MM/*.xyz を渡すこと）。
    out_dir : str
        UMA 最適化後の xyz を書き出すフォルダ。
    cfg :
        UMAConfig 相当のオブジェクト。
        必須フィールド:
            - model_name または model_tag（どちらかあればOK）
            - device : "cpu" or "cuda"
            - fmax : float   （収束閾値 [eV/Å], max|F| 基準）
            - max_steps : int
        任意フィールド:
            - task_name : str  （無ければ "omol"）
            - debug : bool     （情報として読むだけ）

    id_dir : str, optional
        ID ごとのフォルダ（logs/ 以下に log を書き出す）。
    is_abort_now : callable, optional
        即時中断フラグを返す関数。True が返ればその時点で処理を打ち切る。
    chunk_steps : int
        互換性のためのダミー引数（実装上はほぼ使用しない）。

    Returns
    -------
    results : list of dict
        各構造ごとの結果辞書。
        少なくとも以下のキーを持つ:
          - idx
          - xyz_in
          - xyz_out
          - status: "ok" | "error" | "aborted"
          - energy_init_eV
          - energy_final_eV
          - rmsd_final_A
          - message
    """
    if is_abort_now is None:

        def is_abort_now() -> bool:
            return False

    os.makedirs(out_dir, exist_ok=True)

    # cfg.model_name / cfg.model_tag のどちらでも動くようにする
    model_name = getattr(
        cfg,
        "model_name",
        getattr(cfg, "model_tag", "uma-s-1p1"),  # どちらもなければデフォルト
    )
    device = getattr(cfg, "device", "cpu")
    task_name = getattr(cfg, "task_name", "omol")
    fmax_thresh = float(getattr(cfg, "fmax", 0.01))
    max_steps = int(getattr(cfg, "max_steps", 500))
    debug_flag = bool(getattr(cfg, "debug", False))  # 情報として読むだけ

    # ログ出力先ディレクトリ（IDごとの logs/）
    if id_dir:
        logs_dir = os.path.join(id_dir, "logs")
    else:
        # id_dir が無い場合は out_dir を fallback とする
        logs_dir = out_dir
    os.makedirs(logs_dir, exist_ok=True)

    batch_start = time.time()

    # Calculator 構築
    calc = _build_calculator(model_name, device, task_name)

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
            f"# model={model_name}, device={device}, fmax={fmax_thresh}, max_steps={max_steps}, debug={debug_flag}",
        )
        _log_uma_line(
            log_path,
            f"# start_time {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(file_start))}",
        )
        _log_uma_line(log_path, "# step  energy(eV)           f_max(eV/A)       f_rms(eV/A)")

        try:
            atoms = read(xyz_path)
            # charge / spin が無い場合のデフォルト
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
    batch_log = os.path.join(logs_dir, "UMA_batch.log")
    _log_uma_line(
        batch_log,
        f"UMA batch end: n_results={len(results)}, wall_time={batch_wall:.3f} s",
    )

    return results
