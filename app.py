from __future__ import annotations

import os
import time
import shutil
from datetime import datetime, timedelta
from dataclasses import replace
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd
import streamlit as st

from conformer_app.core import (
    AppConfig,
    StageState,
    JobProgressRow,
)
from conformer_app.mm_layer import run_mm_for_job
from conformer_app.uma_layer import run_uma_for_job
from conformer_app.qm_layer import run_qm_for_job
from conformer_app.app_io.excel_io import (
    parse_requirement_excel,
    write_progress_excel,
    create_initial_progress_excel,
    load_progress_excel,
)
from conformer_app.ui import (
    ui_embed_block,
    ui_mm_block,
    ui_uma_block,      # ← UMA 設定（ON/OFF 含む）
    ui_gamess_block,
    ui_output_block,
)
from conformer_app.rdkit.fragments_logic import run_fragments_for_excel
from conformer_app.rdkit.descriptors_logic import run_descriptors_for_excel
from conformer_app.extract.extract import run_geom_qm_extract
from conformer_app.rdkit.db_update_logic import update_batches_with_geom


# =========================
# フラグ・パス設定
# =========================

ABORT_FLAG = Path("./abort_now.flag")
OUTPUT_BASE = Path("./output")
RESULT_ROOT = Path("./result")
TIMING_BASE = RESULT_ROOT

# DB 更新用ファイルの最終出力先（必要に応じてここを書き換える）
DB_EXPORT_ROOT = Path(
    r"\\ad.fujimiinc.co.jp\FUJIMI\CMP事業本部\非公開\CMP開発部_CMP開発課"
    r"\004-実験計画・実験データ\01-実験計画書・関連データ\99-基礎研究"
    r"\R014_MI\93_MI-CS\計算結果_更新用"
)

TOTAL_TIME_LOG = Path("./total_run_time.log")


# =========================
# 即時中断フラグ
# =========================

def abort_now() -> bool:
    """即時中断フラグが立っているかどうか"""
    return ABORT_FLAG.exists()


def clear_abort_flag() -> None:
    """計算終了後などにフラグファイルを削除する"""
    try:
        if ABORT_FLAG.exists():
            ABORT_FLAG.unlink()
    except Exception:
        pass


# =========================
# Formal_Charge の取得（progress.xlsx から）
# =========================

def get_formal_charge_from_progress(
    mol_id: str,
    progress_path: str | Path = "./progress.xlsx",
) -> int | None:
    """
    progress.xlsx の Formal_Charge 列から、ID の形式電荷を取得する。

    - 取得できれば int
    - 失敗・未設定なら None
    """
    path = Path(progress_path)
    if not path.is_file():
        return None

    try:
        df = load_progress_excel(path)
    except Exception:
        return None

    if "ID" not in df.columns or "Formal_Charge" not in df.columns:
        return None

    mol_id_str = str(mol_id).strip()
    mask = df["ID"].astype(str).str.strip() == mol_id_str
    if not mask.any():
        return None

    try:
        val = df.loc[mask, "Formal_Charge"].iloc[0]
        if pd.isna(val):
            return None
        return int(val)
    except Exception:
        return None


# =========================
# QM ポスト処理対象フォルダ探索
# =========================

def find_qm_target_dir(mol_id: str) -> Optional[Path]:
    """
    output/<ID>/QM/gms* のうち、
    「.out と *_final.xyz が両方ある gms フォルダ」の中から
    最も新しいものを 1 つ返す。
    見つからなければ None。
    """
    qm_root = OUTPUT_BASE / str(mol_id) / "QM"
    qm_root = qm_root.resolve()
    if not qm_root.is_dir():
        return None

    candidates: list[Path] = []
    for sub in qm_root.glob("gms*"):
        if not sub.is_dir():
            continue
        outs = list(sub.glob("*.out"))
        xyzs = list(sub.glob("*_final.xyz"))
        if outs and xyzs:
            candidates.append(sub)

    if not candidates:
        return None

    # 更新日時が新しいフォルダを優先
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


# =========================
# タイミングログ
# =========================

def _fmt_sec(dt: Optional[float]) -> str:
    if dt is None:
        return "None"
    return f"{dt:.3f}"


def _fmt_hms_and_sec(dt: Optional[float]) -> str:
    """
    秒数 dt から "HH:MM:SS (xxxx.xxx sec)" という表示文字列を返す。
    """
    if dt is None:
        return "None"
    td = timedelta(seconds=float(dt))
    return f"{td} ({dt:.3f} sec)"


def write_id_timing_log(
    mol_id: str,
    t_start: float,
    t_mm_start: Optional[float],
    t_mm_end: Optional[float],
    t_uma_start: Optional[float],
    t_uma_end: Optional[float],
    t_qm_start: Optional[float],
    t_qm_end: Optional[float],
    note: Optional[str] = None,
) -> None:
    """
    各 ID についてのタイミング情報を
    result/<ID>/run_time.log に追記する。
    """
    log_dir = TIMING_BASE / str(mol_id)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "run_time.log"

    mm_dt = (t_mm_end - t_mm_start) if (t_mm_start is not None and t_mm_end is not None) else None
    uma_dt = (t_uma_end - t_uma_start) if (t_uma_start is not None and t_uma_end is not None) else None
    qm_dt = (t_qm_end - t_qm_start) if (t_qm_start is not None and t_qm_end is not None) else None

    end_candidates = [t for t in (t_qm_end, t_uma_end, t_mm_end) if t is not None]
    total_dt = None
    if end_candidates:
        total_dt = max(end_candidates) - t_start

    start_str = datetime.fromtimestamp(t_start).strftime("%Y-%m-%d %H:%M:%S")

    with open(log_path, "a", encoding="utf-8") as f:
        f.write(f"=== ID {mol_id} ===\n")
        f.write(f"Start_time = {start_str}\n")
        f.write(f"total_time = {_fmt_hms_and_sec(total_dt)}\n")
        f.write(f"MM_time    = {_fmt_hms_and_sec(mm_dt)}\n")
        f.write(f"UMA_time   = {_fmt_hms_and_sec(uma_dt)}\n")
        f.write(f"QM_time    = {_fmt_hms_and_sec(qm_dt)}\n")
        if note:
            f.write(f"note       = {note}\n")
        f.write("\n")


def write_total_time_log(
    total_start: float,
    total_end: float,
    n_jobs: int,
) -> None:
    """
    全 ID の計算にかかった総時間をルート直下 total_run_time.log に追記。
    """
    elapsed = total_end - total_start
    start_str = datetime.fromtimestamp(total_start).strftime("%Y-%m-%d %H:%M:%S")
    end_str = datetime.fromtimestamp(total_end).strftime("%Y-%m-%d %H:%M:%S")

    line = (
        f"Start: {start_str}  "
        f"End: {end_str}  "
        f"Elapsed: {_fmt_hms_and_sec(elapsed)}  "
        f"Jobs: {n_jobs}\n"

    )

    with TOTAL_TIME_LOG.open("a", encoding="utf-8") as f:
        f.write(line)


# =========================
# メインアプリ
# =========================

def main():
    # OpenMP の二重ロード回避（fairchem / torch 用）
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    # セッション状態：プログレス値と「実行中フラグ」を保持
    if "progress_value" not in st.session_state:
        st.session_state["progress_value"] = 0.0
    if "is_running" not in st.session_state:
        st.session_state["is_running"] = False

    st.set_page_config(page_title="コンフォマー自動計算アプリ", layout="wide")
    st.title("DB更新用 自動計算アプリ")
    st.title("～MM / UMA / QM 自動計算 + DBバッチ更新～")

    # ===== サイドバー：進捗 / 中断 =====
    st.sidebar.header("進捗 / 中断")

    if st.sidebar.button("即時中断", type="secondary"):
        ABORT_FLAG.write_text("1", encoding="utf-8")

    if abort_now():
        st.sidebar.warning("即時中断フラグ: ON（UMA / ID 間でチェックされます）")
    else:
        st.sidebar.info("即時中断フラグ: OFF")

    progress_bar = st.sidebar.progress(st.session_state.get("progress_value", 0.0))
    sidebar_status = st.sidebar.empty()

    # ===== 1. Excel 入力 =====
    st.markdown("### 1. Excel 入力")
    st.caption("requirement.xlsx（新規） または progress.xlsx（再開）を指定")

    uploaded = st.file_uploader("requirement / progress Excel (.xlsx)", type=["xlsx"])

    ids: List[str] = []
    smis: List[str] = []
    resume_mode: bool = False

    if uploaded is not None:
        # まず DataFrame として中身を覗いて、requirement か progress か判定
        try:
            df_up = pd.read_excel(uploaded)
        except Exception as e:
            st.error(f"Excel の読み込みに失敗しました: {e}")
            df_up = None

        if df_up is not None:
            has_status_cols = {
                "MM_status",
                "UMA_status",
                "QM_status",
            }.issubset(df_up.columns)

            if has_status_cols:
                # ======================
                # 既存 progress.xlsx (resume モード)
                # ======================
                resume_mode = True

                # アップロードされたファイルをそのまま progress.xlsx として保存
                df_up.to_excel("./progress.xlsx", index=False)

                # ID / SMILES を取得（ヘッダー行 "SMILES" はスキップ）
                col_id = df_up.iloc[:, 0].astype(str).str.strip()
                col_smi = df_up.iloc[:, 1].astype(str).str.strip()

                ids = []
                smis = []
                for id_raw, smi_raw in zip(col_id, col_smi):
                    if not smi_raw:
                        continue
                    if smi_raw.lower() == "smiles":
                        continue
                    if not id_raw:
                        id_raw = smi_raw
                    ids.append(id_raw)
                    smis.append(smi_raw)

                st.success(f"読み込んだ分子数: {len(ids)} 件（resume モード）")
                st.info("既存の progress.xlsx を読み込んだため、DONE 行はスキップして再開します。")

            else:
                # ======================
                # requirement.xlsx モード（新規計算）
                # ======================
                uploaded.seek(0)
                try:
                    ids, smis = parse_requirement_excel(uploaded)
                except Exception as e:
                    st.error(f"Excel の読み込みに失敗しました: {e}")
                    ids, smis = [], []
                else:
                    st.success(f"読み込んだ分子数: {len(ids)} 件")

                if ids:
                    # 初期 progress.xlsx を作成
                    create_initial_progress_excel(ids, smis, "./progress.xlsx")
                    st.info("progress.xlsx を生成しました。")

                    # 官能基フラグメント集計
                    try:
                        progress_path = Path("./progress.xlsx").resolve()
                        result_dir = RESULT_ROOT.resolve()
                        out_frag = run_fragments_for_excel(progress_path, result_dir)
                        st.info(f"官能基フラグメントを出力しました: {out_frag}")
                    except Exception as e:
                        st.warning(f"官能基フラグメントの集計に失敗しました: {e}")

                    # RDKit 記述子集計（原本 + DB 更新用バッチ生成の元データ）
                    try:
                        progress_path = Path("./progress.xlsx").resolve()
                        result_dir = RESULT_ROOT.resolve()
                        out_desc = run_descriptors_for_excel(progress_path, result_dir)
                        st.info(f"RDKit 記述子を出力しました: {out_desc}")
                    except Exception as e:
                        st.warning(f"RDKit 記述子の集計に失敗しました: {e}")

    st.markdown("---")
    st.markdown("### 2. 設定")

    # AppConfig デフォルト作成
    app_cfg = AppConfig()

    # 分子生成：最大コンフォマー数 1000 に初期設定
    try:
        app_cfg.embed.num_confs = 1000
    except Exception:
        pass

    # MM → UMA に渡す構造の上限：最大選抜件数 10
    try:
        app_cfg.mm.max_keep = 10
    except Exception:
        pass

    # UI ブロック（mol_gen / MM / UMA / QM / 出力）
    embed_cfg = ui_embed_block(app_cfg.embed)
    mm_cfg = ui_mm_block(app_cfg.mm)
    uma_cfg, use_uma = ui_uma_block(app_cfg.uma)   # ← UMA 設定 + 使用フラグ
    qm_cfg = ui_gamess_block(app_cfg.qm)
    out_cfg = ui_output_block(app_cfg.output)

    app_cfg.embed = embed_cfg
    app_cfg.mm = mm_cfg
    app_cfg.uma = uma_cfg
    app_cfg.qm = qm_cfg
    app_cfg.output = out_cfg

    st.markdown("---")
    st.markdown("### 3. 実行")

    col_run, col_info = st.columns([1, 3])
    with col_run:
        run_button = st.button("計算実行")
    with col_info:
        if use_uma:
            st.caption("※ MM → UMA → QM を ID ごとに直列で処理します。")
        else:
            st.caption("※ MM → QM（UMA スキップ）を ID ごとに直列で処理します。")
            st.caption("※ UMA OFF 時は、post_MM から最大 3 構造のみ QM に回します。")

    # 入力プレビュー
    if ids:
        with st.expander("入力プレビュー (ID / SMILES)", expanded=False):
            df_preview = pd.DataFrame({"ID": ids, "SMILES": smis})
            st.dataframe(df_preview)

    result_placeholder = st.empty()

    # 実行ボタンが押されたら本体処理
    if run_button and not st.session_state["is_running"]:
        if not ids:
            st.error("先に Excel ファイルを読み込んでください。")
            return

        # -------------------------
        # progress.xlsx から DONE ID を読み取り（resume 用）
        # -------------------------
        done_ids: set[str] = set()
        prog_path = Path("./progress.xlsx")
        if prog_path.is_file():
            try:
                df_prog = load_progress_excel(prog_path)
                has_cols = {"ID", "MM_status", "UMA_status", "QM_status"}.issubset(df_prog.columns)
                if has_cols:
                    for _, r in df_prog.iterrows():
                        rid = str(r["ID"]).strip()
                        mm_s = str(r.get("MM_status", "")).strip().upper()
                        uma_s = str(r.get("UMA_status", "")).strip().upper()
                        qm_s = str(r.get("QM_status", "")).strip().upper()
                        # 3 つとも DONE のものだけ「完全終了」とみなしてスキップ
                        if mm_s == "DONE" and uma_s == "DONE" and qm_s == "DONE":
                            done_ids.add(rid)
            except Exception as e:
                st.warning(f"progress.xlsx の読み取りに失敗したため、resume 無効: {e}")

        # 実際に今回まわす ID/SMILES のペア
        id_smi_pairs: List[Tuple[str, str]] = []
        for mol_id, smi in zip(ids, smis):
            key = str(mol_id).strip()
            if key in done_ids:
                continue
            id_smi_pairs.append((mol_id, smi))

        if not id_smi_pairs:
            st.info("全ての ID が既に DONE のため、新たに実行するジョブはありません。")
            return

        st.session_state["is_running"] = True
        st.session_state["progress_value"] = 0.0
        clear_abort_flag()

        rows: List[JobProgressRow] = []
        n_total = len(id_smi_pairs)

        # 全体実行開始時刻
        total_start = time.time()

        try:
            for idx, (mol_id, smi) in enumerate(id_smi_pairs, start=1):
                row = JobProgressRow(id=mol_id, smiles=smi)

                # タイミング計測用
                t_start = time.time()
                t_mm_start = t_mm_end = None
                t_uma_start = t_uma_end = None
                t_qm_start = t_qm_end = None
                note_for_log: Optional[str] = None

                # Formal charge
                formal_charge = get_formal_charge_from_progress(mol_id)
                if formal_charge is not None:
                    sidebar_status.write(f"[{mol_id}] Formal_Charge = {formal_charge} を使用します。")
                else:
                    sidebar_status.write(f"[{mol_id}] Formal_Charge が見つからないため、UI の電荷設定を使用します。")

                try:
                    # 開始時点で abort されていたらスキップ
                    if abort_now():
                        row.mm_status = StageState.ERROR
                        row.uma_status = StageState.ERROR
                        row.qm_status = StageState.ERROR
                        row.error_message = "即時中断フラグが立っていたため、この ID はスキップされました。"
                        note_for_log = row.error_message
                        write_id_timing_log(
                            mol_id,
                            t_start,
                            t_mm_start, t_mm_end,
                            t_uma_start, t_uma_end,
                            t_qm_start, t_qm_end,
                            note=note_for_log,
                        )
                        rows.append(row)
                        break

                    # ===== MM =====
                    sidebar_status.write(f"[{mol_id}] MM 最適化中...")
                    row.mm_status = StageState.RUN
                    t_mm_start = time.time()
                    mm_result = run_mm_for_job(mol_id, smi, app_cfg.embed, app_cfg.mm, app_cfg.output)
                    t_mm_end = time.time()
                    row.mm_status = StageState.DONE

                    # 進捗テーブル用
                    row.generated_confs = getattr(mm_result, "n_generated", None)
                    row.mm_optimized_confs = getattr(mm_result, "n_optimized", None)
                    row.mm_candidates = getattr(mm_result, "n_candidates", None)
                    row.mm_selected = getattr(mm_result, "n_selected", None)

                    if abort_now():
                        row.uma_status = StageState.ERROR
                        row.qm_status = StageState.ERROR
                        row.error_message = "MM 完了後に中断フラグが立ったため、この ID で停止しました。"
                        note_for_log = row.error_message
                        write_id_timing_log(
                            mol_id,
                            t_start,
                            t_mm_start, t_mm_end,
                            t_uma_start, t_uma_end,
                            t_qm_start, t_qm_end,
                            note=note_for_log,
                        )
                        rows.append(row)
                        write_progress_excel(rows, "./progress.xlsx")
                        break

                    # ===== UMA =====
                    uma_summary_for_qm: Optional[object] = None

                    if use_uma:
                        sidebar_status.write(f"[{mol_id}] UMA 最適化中...")

                        if formal_charge is not None:
                            uma_cfg_for_id = replace(app_cfg.uma, charge=formal_charge)
                        else:
                            uma_cfg_for_id = app_cfg.uma

                        row.uma_status = StageState.RUN
                        t_uma_start = time.time()
                        uma_result = run_uma_for_job(
                            mol_id=mol_id,
                            mm_result=mm_result,       # ★ MM の結果から UMA 入力を決定
                            uma_cfg=uma_cfg_for_id,
                            out_cfg=app_cfg.output,
                            is_abort_now=abort_now,
                        )
                        t_uma_end = time.time()
                        row.uma_status = StageState.DONE

                        # ★ UMAResult をそのまま QM に渡す（UMA ありルートを確実に使う）
                        uma_summary_for_qm = uma_result

                        if abort_now():
                            row.qm_status = StageState.ERROR
                            row.error_message = "UMA 完了後に中断フラグが立ったため、QM をスキップして停止しました。"
                            note_for_log = row.error_message
                            write_id_timing_log(
                                mol_id,
                                t_start,
                                t_mm_start, t_mm_end,
                                t_uma_start, t_uma_end,
                                t_qm_start, t_qm_end,
                                note=note_for_log,
                            )
                            rows.append(row)
                            write_progress_excel(rows, "./progress.xlsx")
                            break
                    else:
                        # UMA を使わない場合は、post_MM から QM 構造を準備するので
                        # row.uma_status は DONE 扱いにしておき、QM 側には None を渡す
                        row.uma_status = StageState.DONE
                        uma_summary_for_qm = None
                        note_for_log = "UMA disabled for this job."

                    # ===== QM =====
                    sidebar_status.write(f"[{mol_id}] GAMESS 実行中...")

                    if formal_charge is not None:
                        qm_cfg_for_id = replace(app_cfg.qm, icharge=formal_charge)
                    else:
                        qm_cfg_for_id = app_cfg.qm

                    row.qm_status = StageState.RUN
                    t_qm_start = time.time()
                    _qm_result = run_qm_for_job(
                        mol_id,
                        uma_summary_for_qm,   # UMA あり: UMAResult / なし: None
                        qm_cfg_for_id,
                        app_cfg.output,
                    )
                    t_qm_end = time.time()
                    row.qm_status = StageState.DONE

                    # ===== Geom + QM 統合ポスト処理 =====
                    try:
                        qm_target_dir = find_qm_target_dir(str(mol_id))
                        if qm_target_dir is None:
                            msg = (
                                f"[{mol_id}] output/{mol_id}/QM/gms* 内に "
                                f".out と *_final.xyz を持つフォルダが見つからなかったため、"
                                f"Geom+QM 統合をスキップしました。"
                            )
                            sidebar_status.warning(msg)
                            note_for_log = msg
                        else:
                            sidebar_status.write(f"[{mol_id}] Geom+QM 統合解析中...")
                            geom_out_path = run_geom_qm_extract(qm_target_dir)

                            # Geom_QM_summary.xlsx を result/<ID> にコピー
                            result_id_dir = RESULT_ROOT / str(mol_id)
                            result_id_dir.mkdir(parents=True, exist_ok=True)

                            if geom_out_path is not None and geom_out_path.is_file():
                                final_path = result_id_dir / geom_out_path.name
                                shutil.copy2(geom_out_path, final_path)
                                msg = f"[{mol_id}] Geom_QM_summary.xlsx を {final_path} に出力しました。"
                                sidebar_status.write(msg)
                                note_for_log = msg
                            else:
                                msg = (
                                    f"[{mol_id}] Geom_QM_summary.xlsx が生成されませんでした。"
                                    f" (.out / *_final.xyz が不足している可能性があります)"
                                )
                                sidebar_status.warning(msg)
                                note_for_log = msg

                    except Exception as e:
                        msg = f"[{mol_id}] Geom+QM 統合ポスト処理でエラー: {e}"
                        sidebar_status.warning(msg)
                        note_for_log = msg

                    # ===== DB 更新用バッチの更新 =====
                    try:
                        update_batches_with_geom(
                            RESULT_ROOT.resolve(),
                            export_root=DB_EXPORT_ROOT.resolve(),
                        )
                    except Exception as e:
                        sidebar_status.warning(
                            f"[{mol_id}] DB 更新用バッチ更新中にエラー: {e}"
                        )

                    sidebar_status.write(f"[{mol_id}] 完了")

                except Exception as e:
                    # ステータスを ERROR に更新
                    if row.mm_status == StageState.RUN:
                        row.mm_status = StageState.ERROR
                    elif row.uma_status == StageState.RUN:
                        row.uma_status = StageState.ERROR
                    elif row.qm_status == StageState.RUN:
                        row.qm_status = StageState.ERROR
                    row.error_message = f"実行中にエラーが発生しました: {e}"
                    sidebar_status.error(f"[{mol_id}] エラー: {e}")
                    note_for_log = row.error_message

                finally:
                    # 各 ID ごとにタイミングログを出力
                    write_id_timing_log(
                        mol_id,
                        t_start,
                        t_mm_start, t_mm_end,
                        t_uma_start, t_uma_end,
                        t_qm_start, t_qm_end,
                        note=note_for_log,
                    )

                # 結果を蓄積・Excel 反映・進捗更新
                rows.append(row)
                write_progress_excel(rows, "./progress.xlsx")
                st.session_state["progress_value"] = idx / n_total
                progress_bar.progress(st.session_state["progress_value"])

        finally:
            clear_abort_flag()
            st.session_state["is_running"] = False

            # 全体時間ログ
            total_end = time.time()
            try:
                write_total_time_log(total_start, total_end, n_total)
            except Exception:
                pass

        # 最終結果を画面表示
        df_res = pd.DataFrame(
            [
                {
                    "ID": r.id,
                    "SMILES": r.smiles,
                    "MM_status": r.mm_status.value,
                    "UMA_status": r.uma_status.value,
                    "QM_status": r.qm_status.value,
                    "error_message": r.error_message,
                }
                for r in rows
            ]
        )
        result_placeholder.dataframe(df_res)

        st.success(
            "計算が終了しました（または中断されました）。"
            "progress.xlsx / output / result / DB_EXPORT_ROOT / total_run_time.log を確認してください。"
        )
        sidebar_status.write("全ジョブ終了")


if __name__ == "__main__":
    main()
