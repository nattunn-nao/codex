from __future__ import annotations

import os
from typing import Any, Dict, List

import pandas as pd


def write_uma_summary(id_dir: str, results: List[Dict[str, Any]]) -> str:
    """
    UMA 最適化結果を ID フォルダ直下に UMA_summary.xlsx として保存する。

    Parameters
    ----------
    id_dir : str
        1つの ID に対応するフォルダ（例: ./output/111-1）。
    results : list of dict
        uma_optimize_xyz_batch_fairchem() が返す結果リスト。

    Returns
    -------
    out_path : str
        出力した UMA_summary.xlsx のパス（失敗した場合は空文字）。
    """
    if not results:
        return ""

    try:
        df = pd.DataFrame(results)
    except Exception:
        return ""

    out_path = os.path.join(id_dir, "UMA_summary.xlsx")
    try:
        df.to_excel(out_path, index=False, sheet_name="UMA")
    except Exception:
        return ""

    return out_path
