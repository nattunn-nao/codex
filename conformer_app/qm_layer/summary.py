# conformer_app/qm_layer/summary.py
from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Dict, List

import pandas as pd


def _log_side(msg: str) -> None:
    try:
        import streamlit as st

        st.sidebar.write(msg)
    except Exception:
        pass


def write_qm_summary(id_dir: str, records: List[Dict[str, object]]) -> str:
    """
    gmsX ごとの結果レコードから <ID>/QM_summary.xlsx を作成する。
    """
    if not records:
        _log_side(f"[{os.path.basename(id_dir)}] QM records が空のため QM_summary.xlsx は出力されません。")
        return ""

    df = pd.DataFrame(records)
    out_path = Path(id_dir) / "QM_summary.xlsx"

    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_excel(out_path, index=False, sheet_name="QM")
        _log_side(f"[{os.path.basename(id_dir)}] QM_summary.xlsx を出力: {out_path}")
    except Exception as e:
        _log_side(f"[{os.path.basename(id_dir)}] QM_summary.xlsx の書き込みに失敗: {e}")
        return ""

    return str(out_path)


def copy_best_structure(
    id_dir: str,
    mol_id: str,
    qm_summary_path: str,
    result_root: str = "./result",
) -> str:
    """
    QM_summary.xlsx を参照し、エネルギー最小の構造を ./result/<ID>/best.xyz にコピーする。

    何かの理由で best が決まらない場合でも、
      - result_root
      - result_root/<ID>
    までは作成しておく。
    """
    basename = os.path.basename(id_dir)

    # result_root / <ID> はとりあえず作成しておく
    result_dir = Path(result_root) / mol_id
    try:
        result_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        _log_side(f"[{basename}] result ディレクトリ作成に失敗: {e}")
        return ""

    if not qm_summary_path:
        _log_side(f"[{basename}] QM_summary.xlsx のパスが空のため best 構造を決定できません。")
        return ""

    if not os.path.exists(qm_summary_path):
        _log_side(f"[{basename}] QM_summary.xlsx が見つかりません: {qm_summary_path}")
        return ""

    try:
        df = pd.read_excel(qm_summary_path)
    except Exception as e:
        _log_side(f"[{basename}] QM_summary.xlsx の読み込みに失敗: {e}")
        return ""

    if "status" not in df.columns or "final_xyz" not in df.columns:
        _log_side(f"[{basename}] QM_summary.xlsx のカラムに status / final_xyz が見つかりません。")
        return ""

    df_use = df[(df["status"] == "ok") & df["final_xyz"].notna()]
    if df_use.empty:
        _log_side(f"[{basename}] status=='ok' かつ final_xyz が有効なレコードがありません。")
        return ""

    # エネルギーでソート（あれば）
    if "energy_hartree" in df_use.columns:
        df_use = df_use.sort_values("energy_hartree", ascending=True)

    row0 = df_use.iloc[0]
    src_xyz = str(row0["final_xyz"])
    if not os.path.exists(src_xyz):
        _log_side(f"[{basename}] best 候補の final_xyz が見つかりません: {src_xyz}")
        return ""

    dst_xyz = result_dir / "best.xyz"

    try:
        shutil.copy2(src_xyz, dst_xyz)
        _log_side(f"[{basename}] 最安定構造を result へコピー: {dst_xyz}")
    except Exception as e:
        _log_side(f"[{basename}] best.xyz コピーに失敗: {e}")
        return ""

    return str(dst_xyz)
