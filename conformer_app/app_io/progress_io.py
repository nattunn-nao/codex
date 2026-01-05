# conformer_app/app_io/progress_io.py
from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, List, Tuple

import pandas as pd


# Progress.xlsx の標準カラム
PROGRESS_COLUMNS = ["ID", "SMILES", "MM_status", "UMA_status", "QM_status", "error_message"]


def _empty_progress_df() -> pd.DataFrame:
    """空の Progress DataFrame を作る。"""
    df = pd.DataFrame(columns=PROGRESS_COLUMNS)
    return df


def _normalize_progress_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    既存の Progress.xlsx を読み込んだときに、
    必要カラムをすべて持つように補正し、カラム順も揃える。
    """
    df = df.copy()

    # 必須カラムの追加
    if "ID" not in df.columns:
        df["ID"] = ""
    if "SMILES" not in df.columns:
        df["SMILES"] = ""
    if "MM_status" not in df.columns:
        df["MM_status"] = "WAIT"
    if "UMA_status" not in df.columns:
        df["UMA_status"] = "WAIT"
    if "QM_status" not in df.columns:
        df["QM_status"] = "WAIT"
    if "error_message" not in df.columns:
        df["error_message"] = ""

    # カラム順を固定
    df = df[PROGRESS_COLUMNS]
    return df


def load_progress(progress_path: str | os.PathLike) -> pd.DataFrame:
    """
    Progress.xlsx を読み込んで標準形に揃えた DataFrame を返す。
    ファイルが無ければ空の Progress を返す。
    """
    path = Path(progress_path)
    if not path.is_file():
        return _empty_progress_df()

    try:
        df = pd.read_excel(path)
    except Exception:
        # 壊れているなどで読み込めなければ空で作り直す
        return _empty_progress_df()

    return _normalize_progress_df(df)


def save_progress(df: pd.DataFrame, progress_path: str | os.PathLike) -> None:
    """
    Progress DataFrame を Excel として保存する。
    """
    path = Path(progress_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df = _normalize_progress_df(df)
    df.to_excel(path, index=False)


def ensure_progress_for_entries(
    entries: Iterable[Tuple[str, str, str]],
    progress_path: str | os.PathLike = "./progress.xlsx",
) -> pd.DataFrame:
    """
    entries で与えられた全 ID / SMILES を Progress.xlsx に反映する。

    Parameters
    ----------
    entries : Iterable[Tuple[str, str, str]]
        (smiles, mol_id, raw_id) のタプル列を想定。
        Progress ファイル上では mol_id を ID、smiles を SMILES として用いる。
        raw_id はここでは使わないが、呼び出し側の互換性のために受け取る。
    progress_path : str or Path
        Progress.xlsx のパス。

    挙動
    ----
    - 既存の Progress.xlsx があれば読み込んでから更新する。
    - entries に含まれる **すべての mol_id** が Progress に 1 行ずつ存在するようにする。
    - 既に同じ ID の行があれば、その行の SMILES を最新の値で上書きする。
    - 新規に追加される行のステータスは
        MM_status / UMA_status / QM_status = "WAIT"
        error_message = ""
      で初期化する。
    - 既に Progress に存在している ID のステータスはそのまま維持する。
    """
    df = load_progress(progress_path)

    # ID -> 行 index のマップを作成（文字列 ID で管理）
    if "ID" in df.columns:
        id_series = df["ID"].astype(str)
        id_to_idx = {id_: i for i, id_ in enumerate(id_series)}
    else:
        id_to_idx = {}

    for smi, mol_id, raw_id in entries:
        mol_id_str = str(mol_id)
        smi_str = str(smi)

        if mol_id_str in id_to_idx:
            # 既存行を更新：SMILES は常に最新の値で上書きする
            i = id_to_idx[mol_id_str]
            df.at[i, "ID"] = mol_id_str
            df.at[i, "SMILES"] = smi_str
        else:
            # 新規行を追加：ステータスはすべて WAIT で初期化
            new_row = {
                "ID": mol_id_str,
                "SMILES": smi_str,
                "MM_status": "WAIT",
                "UMA_status": "WAIT",
                "QM_status": "WAIT",
                "error_message": "",
            }
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
            id_to_idx[mol_id_str] = len(df) - 1

    df = _normalize_progress_df(df)
    save_progress(df, progress_path)
    return df
