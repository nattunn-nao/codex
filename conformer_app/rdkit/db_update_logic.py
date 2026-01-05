from __future__ import annotations

"""
conformer_app/rdkit/db_update_logic.py

- result/progress_rdkit_descriptors.xlsx を原本として読み込み
- ID ごとに result/<ID>/Geom_QM_summary.xlsx を探して幾何・QM情報をマージ
- batch_size 個ずつ `db_update_work/db_batch_XXX.xlsx` を更新
- そのバッチ内の全 ID に Geom_QM_summary がそろったら
  `db_update_export/db_batch_XXX.xlsx` にコピー（DB 取り込み用）

追加で算出する列（単位 eV）:
    dE_HOMO_LUMO   = LUMO_eV      - HOMO_eV
    dE_next_HOMO   = HOMO_eV      - next_HOMO_eV
    dE_next_LUMO   = next_LUMO_eV - LUMO_eV

SASA 関連:
    polar_SASA_ratio = SASA_polar_area_A2 / SASA_area_A2
"""

from pathlib import Path
from typing import List, Dict, Any, Optional

import pandas as pd

DESCRIPTOR_FILENAME = "progress_rdkit_descriptors.xlsx"
WORK_DIR_NAME = "db_update_work"
EXPORT_DIR_NAME = "db_update_export"

BATCH_SIZE_DEFAULT = 10

# Geom_QM_summary から直接コピーする列
GEOM_BASE_COLS = [
    "dipole_moment",
    "HOMO_eV",
    "LUMO_eV",
    "next_HOMO_eV",
    "next_LUMO_eV",
    "SASA_area_A2",
    "SASA_polar_area_A2",
    "vdW_area_A2",
    "vdW_volume_A3",
    "Rg_A",
    "ovality",
]

# ここで新しく計算する列
GEOM_DERIVED_COLS = [
    "dE_HOMO_LUMO",
    "dE_next_HOMO",
    "dE_next_LUMO",
    "polar_SASA_ratio",
]

ALL_GEOM_COLS = GEOM_BASE_COLS + GEOM_DERIVED_COLS


# ============================================================
# 内部ユーティリティ
# ============================================================


def _load_descriptors(result_root: Path) -> pd.DataFrame:
    """
    result_root/progress_rdkit_descriptors.xlsx を読み込む。
    ID は文字列として扱う。
    """
    desc_path = result_root / DESCRIPTOR_FILENAME
    if not desc_path.is_file():
        raise FileNotFoundError(f"記述子 Excel が見つかりません: {desc_path}")

    df = pd.read_excel(desc_path)
    if "ID" not in df.columns:
        raise ValueError(f"'ID' 列が見つかりません: {desc_path}")

    df["ID"] = df["ID"].astype(str).str.strip()
    return df


def _ensure_geom_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    df に幾何・QM 用の列がなければ追加する（初期値 None）。
    """
    for col in ALL_GEOM_COLS:
        if col not in df.columns:
            df[col] = None
    return df


def _load_geom_row(result_root: Path, mol_id: str) -> Optional[pd.Series]:
    """
    result_root/<ID>/Geom_QM_summary.xlsx の 1 行目を返す。
    なければ None。
    """
    geom_path = result_root / str(mol_id) / "Geom_QM_summary.xlsx"
    if not geom_path.is_file():
        return None

    try:
        gdf = pd.read_excel(geom_path)
    except Exception:
        return None

    if gdf.empty:
        return None

    # 通常 1 行だけのはずなので先頭行を返す
    return gdf.iloc[0]


def _safe_delta(a: Any, b: Any) -> Optional[float]:
    """
    b - a を安全に計算（どちらか NaN/None なら None）。
    """
    try:
        if pd.isna(a) or pd.isna(b):
            return None
        return float(b) - float(a)
    except Exception:
        return None


def _fill_geom_for_batch(df_batch: pd.DataFrame, result_root: Path) -> pd.DataFrame:
    """
    1 バッチぶんの DataFrame に対して、
    result_root/<ID>/Geom_QM_summary.xlsx を見て幾何・QM 情報を埋める。
    """
    df_batch = _ensure_geom_columns(df_batch)

    for idx, row in df_batch.iterrows():
        mol_id = str(row["ID"]).strip()
        geom_row = _load_geom_row(result_root, mol_id)
        if geom_row is None:
            # まだ Geom_QM_summary が出来ていない ID
            continue

        # --- 直接コピーする列 ---
        for col in GEOM_BASE_COLS:
            if col in geom_row.index:
                df_batch.at[idx, col] = geom_row.get(col)

        # --- エネルギー差の計算（eV） ---
        HOMO = geom_row.get("HOMO_eV")
        next_HOMO = geom_row.get("next_HOMO_eV")
        LUMO = geom_row.get("LUMO_eV")
        next_LUMO = geom_row.get("next_LUMO_eV")

        df_batch.at[idx, "dE_HOMO_LUMO"] = _safe_delta(HOMO, LUMO)
        df_batch.at[idx, "dE_next_HOMO"] = _safe_delta(next_HOMO, HOMO)
        df_batch.at[idx, "dE_next_LUMO"] = _safe_delta(LUMO, next_LUMO)

        # --- polar_SASA_ratio の計算 ---
        SASA = geom_row.get("SASA_area_A2")
        SASA_pol = geom_row.get("SASA_polar_area_A2")
        try:
            if (
                (SASA is not None)
                and (not pd.isna(SASA))
                and float(SASA) != 0.0
                and (SASA_pol is not None)
                and (not pd.isna(SASA_pol))
            ):
                df_batch.at[idx, "polar_SASA_ratio"] = float(SASA_pol) / float(SASA)
            else:
                df_batch.at[idx, "polar_SASA_ratio"] = None
        except Exception:
            df_batch.at[idx, "polar_SASA_ratio"] = None

    return df_batch


def _batch_is_complete(df_batch: pd.DataFrame) -> bool:
    """
    このバッチの全 ID について Geom_QM_summary がそろっているか？
    ここでは「HOMO_eV が全部埋まっているか」で判定する。
    """
    if "HOMO_eV" not in df_batch.columns:
        return False
    return df_batch["HOMO_eV"].notna().all()


def _iter_descriptor_batches(df_desc: pd.DataFrame, batch_size: int):
    """
    記述子 DataFrame を batch_size ごとのチャンクに分割して返すジェネレータ。

    yield: (batch_index (1-origin), df_batch)
    """
    n = len(df_desc)
    if n == 0:
        return

    # batch_size 以下しかなくても 1 つのバッチとして返す
    if batch_size <= 0:
        batch_size = n

    batch_idx = 1
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        yield batch_idx, df_desc.iloc[start:end].copy()
        batch_idx += 1


# ============================================================
# 外部公開関数
# ============================================================


def update_batches_with_geom(
    result_root: Path,
    batch_size: int = BATCH_SIZE_DEFAULT,
    export_root: Optional[Path] = None,
) -> List[Path]:
    """
    DB 更新用バッチを「幾何・QM 結果を含めて」更新する。

    - result_root/progress_rdkit_descriptors.xlsx を読み込み
    - ID を batch_size 個ずつに分割して
        result_root/db_update_work/db_batch_XXX.xlsx
      を作成/更新
    - 各バッチについて、全行の HOMO_eV が埋まっていれば
        result_root/db_update_export/db_batch_XXX.xlsx
      にコピー（export_root を指定した場合はそちら配下）

    Returns
    -------
    List[Path]
        今回「export 側」に書き出したバッチファイルのパス一覧。
    """
    result_root = result_root.resolve()
    if export_root is None:
        export_root = result_root / EXPORT_DIR_NAME
    export_root = export_root.resolve()

    work_root = result_root / WORK_DIR_NAME
    work_root.mkdir(parents=True, exist_ok=True)
    export_root.mkdir(parents=True, exist_ok=True)

    # 1. 記述子の原本を読み込み
    df_desc = _load_descriptors(result_root)

    exported: List[Path] = []

    # 2. バッチごとに Geom_QM_summary をマージして保存
    for batch_idx, df_batch in _iter_descriptor_batches(df_desc, batch_size):
        df_batch = _ensure_geom_columns(df_batch)
        df_batch = _fill_geom_for_batch(df_batch, result_root)

        work_path = work_root / f"db_batch_{batch_idx:03d}.xlsx"
        df_batch.to_excel(work_path, index=False)

        # 3. 全行の HOMO_eV が埋まっていれば export 側にもコピー
        if _batch_is_complete(df_batch):
            export_path = export_root / work_path.name
            df_batch.to_excel(export_path, index=False)
            exported.append(export_path)

    return exported
