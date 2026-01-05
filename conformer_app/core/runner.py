from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import pandas as pd
import streamlit as st

from conformer_app.core.config import (
    AppConfig,
    MMConfig,
    UMAConfig,
    QMConfig,
    OutputConfig,
)
from conformer_app.mm_layer import MMResult, run_mm_for_job
from conformer_app.uma_layer import run_uma_for_job
from conformer_app.qm_layer import QMResult, run_qm_for_job
from conformer_app.app_io.progress_io import ensure_progress_for_entries


# =========================
# 進捗テーブル用の行定義
# =========================

@dataclass
class ProgressRow:
    """progress.xlsx と画面テーブル用の 1 行."""

    ID: str
    SMILES: str

    # カウント系（MM）
    target_confs: int = 0
    generated_confs: int = 0
    mm_optimized_confs: int = 0
    mm_candidates: int = 0
    mm_selected: int = 0

    # ステータス
    MM_status: str = "WAIT"   # WAIT / RUN / DONE / ERROR / ABORTED
    UMA_status: str = "WAIT"  # WAIT / RUN / DONE / ERROR / ABORTED
    QM_status: str = "WAIT"   # WAIT / RUN / PREPARED / DONE / ERROR / ABORTED

    # エラーメッセージ
    error_message: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# =========================
# 内部ユーティリティ
# =========================

def _log_side(msg: str) -> None:
    """Streamlit サイドバーに軽くログ出力。"""
    try:
        st.sidebar.write(msg)
    except Exception:
        pass


def _progress_path_from_out_cfg(out_cfg: OutputConfig) -> str:
    """
    OutputConfig から progress.xlsx のパスを決める。
    なければ ./progress.xlsx にする。
    """
    path = getattr(out_cfg, "progress_path", "") or "progress.xlsx"
    return path


def _write_progress_excel(
    rows: List[ProgressRow],
    out_cfg: OutputConfig,
) -> None:
    """
    rows を progress.xlsx として保存する。
    """
    if not rows:
        return

    df = pd.DataFrame([r.to_dict() for r in rows])
    progress_path = _progress_path_from_out_cfg(out_cfg)
    try:
        Path(progress_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_excel(progress_path, index=False, sheet_name="progress")
    except Exception as e:
        _log_side(f"progress.xlsx の書き込みに失敗: {e}")


def _init_progress_rows_from_progress(
    entries: List[Tuple[str, str, str]],
    out_cfg: OutputConfig,
) -> List[ProgressRow]:
    """
    entries と out_cfg から Progress.xlsx を用意し、それを ProgressRow に変換する。

    - 既存の Progress.xlsx があればそれを読み込む。
    - なければ entries から新しく Progress.xlsx を作成する。
    """
    progress_path = _progress_path_from_out_cfg(out_cfg)
    df = ensure_progress_for_entries(entries, progress_path=progress_path)

    rows: List[ProgressRow] = []
    for _, row in df.iterrows():
        rows.append(
            ProgressRow(
                ID=str(row.get("ID", "")).strip(),
                SMILES=str(row.get("SMILES", "")).strip(),
                MM_status=str(row.get("MM_status", "WAIT")),
                UMA_status=str(row.get("UMA_status", "WAIT")),
                QM_status=str(row.get("QM_status", "WAIT")),
                error_message=str(row.get("error_message", "")),
            )
        )
    return rows


def _find_row(rows: List[ProgressRow], mol_id: str) -> ProgressRow:
    """
    mol_id に対応する ProgressRow を返す。無ければ新規作成して返す。
    """
    for r in rows:
        if r.ID == mol_id:
            return r
    r = ProgressRow(ID=mol_id, SMILES="")
    rows.append(r)
    return r


def _reset_row_for_rerun(row: ProgressRow) -> None:
    """
    ある ID を「頭からやり直す」ときに、その ID の ProgressRow をリセットする。

    - ステータスを WAIT に戻す
    - エラー文字列を消す
    - カウント系も 0 に戻す
    """
    row.MM_status = "WAIT"
    row.UMA_status = "WAIT"
    row.QM_status = "WAIT"
    row.error_message = ""

    row.target_confs = 0
    row.generated_confs = 0
    row.mm_optimized_confs = 0
    row.mm_candidates = 0
    row.mm_selected = 0


# =========================
# メイン: すべての ID を回す
# =========================

def run_all_ids(
    entries: List[Tuple[str, str, str]],
    app_cfg: AppConfig,
    mm_cfg: MMConfig,
    uma_cfg: UMAConfig,
    qm_cfg: QMConfig,
    out_cfg: OutputConfig,
    # 将来用の互換引数（現状は未使用）
    resume_file: Optional[Any] = None,
    # 中断フラグ:
    #   is_abort_now : 「即時中断」ボタン用 → True なら今のステップが終わったら ID ループを即終了
    #   is_wait_stop : 「待機終了」ボタン用 → True なら「現在の ID を最後までやってから」終了
    is_abort_now: Optional[Callable[[], bool]] = None,
    is_wait_stop: Optional[Callable[[], bool]] = None,
) -> None:
    """
    requirement Excel から読み込んだ entries を、ID ごとに
      MM → UMA → QM
    の順番で直列実行するドライバ。

    再開ポリシー:
      - Progress.xlsx を進捗の唯一のソースとする。
      - QM_status == "DONE" の ID は完全にスキップ。
      - それ以外の ID は「MM から UMA, QM まで」を頭からやり直す。
    """
    if is_abort_now is None:
        is_abort_now = lambda: False
    if is_wait_stop is None:
        is_wait_stop = lambda: False

    # 進捗行の初期化（Progress.xlsx を用意 or 読み込み）
    rows: List[ProgressRow] = _init_progress_rows_from_progress(entries, out_cfg)

    total_ids = len(entries)
    finished_ids = 0

    # メインループ: ID ごとに直列
    for smi, mol_id, raw_id in entries:
        row = _find_row(rows, mol_id)

        # すでに QM DONE の ID は完全スキップ
        if row.QM_status == "DONE":
            _log_side(f"[{mol_id}] QM DONE (progress). この ID はスキップします。")
            finished_ids += 1
            _write_progress_excel(rows, out_cfg)
            continue

        # ID 開始前の即時中断チェック
        if is_abort_now():
            _log_side(f"[{mol_id}] 即時中断フラグ検出。処理を中断します。")
            if row.MM_status == "WAIT":
                row.MM_status = "ABORTED"
            if row.UMA_status == "WAIT":
                row.UMA_status = "ABORTED"
            if row.QM_status == "WAIT":
                row.QM_status = "ABORTED"
            row.error_message = "aborted before start"
            _write_progress_excel(rows, out_cfg)
            break

        # この ID は「頭からやり直す」のでステータスとカウントをリセット
        _reset_row_for_rerun(row)
        _write_progress_excel(rows, out_cfg)

        # ========= MM ステージ =========
        mm_result: Optional[MMResult] = None
        try:
            row.MM_status = "RUN"
            _write_progress_excel(rows, out_cfg)

            mm_result = run_mm_for_job(
                mol_id=mol_id,
                smiles=smi,
                app_cfg=app_cfg,
                mm_cfg=mm_cfg,
                out_cfg=out_cfg,
                is_abort_now=is_abort_now,
            )

            # MMResult からカウントを反映
            row.target_confs = getattr(mm_result, "target_confs", row.target_confs)
            row.generated_confs = getattr(mm_result, "generated_confs", row.generated_confs)
            row.mm_optimized_confs = getattr(mm_result, "mm_optimized_confs", row.mm_optimized_confs)
            row.mm_candidates = getattr(mm_result, "mm_candidates", row.mm_candidates)
            row.mm_selected = getattr(mm_result, "mm_selected", row.mm_selected)

            if mm_result.status == "ok":
                row.MM_status = "DONE"
            elif mm_result.status == "aborted":
                row.MM_status = "ABORTED"
                row.error_message = mm_result.message or "MM aborted"
                _write_progress_excel(rows, out_cfg)
                break  # ID ループを抜ける
            else:
                row.MM_status = "ERROR"
                row.error_message = mm_result.message or "MM error"
                _write_progress_excel(rows, out_cfg)
                finished_ids += 1
                continue  # この ID はこれ以上進めない

        except Exception as e:
            row.MM_status = "ERROR"
            row.error_message = f"MM exception: {e}"
            _log_side(f"[{mol_id}] MM で例外発生: {e}")
            _write_progress_excel(rows, out_cfg)
            finished_ids += 1
            continue

        _write_progress_excel(rows, out_cfg)

        # ========= UMA ステージ =========
        uma_enabled = getattr(uma_cfg, "enabled", True)
        uma_result: Any = None

        if uma_enabled:
            try:
                row.UMA_status = "RUN"
                _write_progress_excel(rows, out_cfg)

                uma_result = run_uma_for_job(
                    mol_id=mol_id,
                    mm_result=mm_result,
                    uma_cfg=uma_cfg,
                    out_cfg=out_cfg,
                    is_abort_now=is_abort_now,
                )

                status = getattr(uma_result, "status", "error")
                msg = getattr(uma_result, "message", "")

                if status == "ok":
                    row.UMA_status = "DONE"
                elif status == "aborted":
                    row.UMA_status = "ABORTED"
                    row.error_message = msg or "UMA aborted"
                    _write_progress_excel(rows, out_cfg)
                    break
                else:
                    row.UMA_status = "ERROR"
                    row.error_message = msg or "UMA error"
                    _write_progress_excel(rows, out_cfg)
                    finished_ids += 1
                    continue

            except Exception as e:
                row.UMA_status = "ERROR"
                row.error_message = f"UMA exception: {e}"
                _log_side(f"[{mol_id}] UMA で例外発生: {e}")
                _write_progress_excel(rows, out_cfg)
                finished_ids += 1
                continue

            _write_progress_excel(rows, out_cfg)
        else:
            row.UMA_status = row.UMA_status or "WAIT"
            uma_result = None

        # ========= QM ステージ =========
        qm_enabled = getattr(qm_cfg, "enabled", True)
        qm_result: Optional[QMResult] = None

        if qm_enabled:
            try:
                row.QM_status = "RUN"
                _write_progress_excel(rows, out_cfg)

                # UMA を実行していれば、その summary.xlsx パスを qm_layer に渡す
                uma_summary_path = ""
                if uma_result is not None:
                    uma_summary_path = getattr(uma_result, "uma_summary_path", "")
                if not uma_summary_path:
                    # UMA をスキップした場合・古い実行から再開の場合のフォールバック
                    uma_summary_path = str(
                        Path(out_cfg.base_dir) / mol_id / "UMA_summary.xlsx"
                    )

                qm_result = run_qm_for_job(
                    mol_id=mol_id,
                    uma_summary_path=uma_summary_path,
                    qm_config=qm_cfg,
                    out_cfg=out_cfg,
                    max_jobs=3,
                )

                if qm_result is not None:
                    if getattr(qm_result, "best_xyz_path", ""):
                        row.QM_status = "DONE"
                    else:
                        # summary まではできたが best 決定できなかったケースも DONE 相当扱い
                        row.QM_status = "DONE"
                else:
                    row.QM_status = "ERROR"
                    row.error_message = "QM_result is None"

            except Exception as e:
                row.QM_status = "ERROR"
                row.error_message = f"QM exception: {e}"
                _log_side(f"[{mol_id}] QM で例外発生: {e}")

            _write_progress_excel(rows, out_cfg)
        else:
            row.QM_status = row.QM_status or "WAIT"

        finished_ids += 1
        _write_progress_excel(rows, out_cfg)

        # ===== 「待機終了」フラグチェック =====
        if is_wait_stop():
            _log_side(f"[{mol_id}] 待機終了フラグ検出。この ID までで停止します。")
            break

        # ===== 「即時中断」フラグも一応ここでもチェック =====
        if is_abort_now():
            _log_side(f"[{mol_id}] 即時中断フラグ検出。この ID までで停止します。")
            break

    # ループ終了後、最終版 progress.xlsx を出力
    _write_progress_excel(rows, out_cfg)
    _log_side("run_all_ids が終了しました。")
