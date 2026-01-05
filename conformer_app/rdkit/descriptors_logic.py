from __future__ import annotations

"""
conformer_app/rdkit/descriptors_logic.py

progress.xlsx（A列: ID, B列: SMILES）から
RDKit 記述子を計算し、

    result/progress_rdkit_descriptors.xlsx

という原本 Excel を出力するモジュール。

- ここでは「IDごとの RDKit 記述子の集計」だけを行う。
- DB 更新用のバッチ分割や Geom_QM_summary.xlsx とのマージは、
  conformer_app.rdkit.db_update_logic.update_batches_with_geom()
  を app.py 側から呼び出して実施する。
"""

from pathlib import Path
from typing import Dict

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, Crippen


# =========================
# 1. Excel 入出力まわり
# =========================


def load_id_smiles_from_excel(path: Path) -> pd.DataFrame:
    """
    Excel から A列(0), B列(1)だけを読み込んで DataFrame を返す。
    1行目はヘッダー、2行目以降がデータを想定。

    読み込んだ列名は強制的に ["ID", "SMILES"] にする。
    """
    df = pd.read_excel(path, usecols=[0, 1])
    df.columns = ["ID", "SMILES"]
    # ID は文字列として扱う
    df["ID"] = df["ID"].astype(str).str.strip()
    return df


def make_output_path(result_root: Path) -> Path:
    """
    result_root（例: ./result）配下に、
    progress_rdkit_descriptors.xlsx のパスを返す。
    """
    result_root.mkdir(parents=True, exist_ok=True)
    return result_root / "progress_rdkit_descriptors.xlsx"


# =========================
# 2. RDKit: 1分子ぶんの記述子計算
# =========================


def smiles_to_mol(smiles: str):
    """SMILES 文字列から RDKit Mol を作る。失敗したら None を返す。"""
    if pd.isna(smiles):
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    return mol


def calc_descriptors_for_mol(mol) -> Dict[str, float | int | None]:
    """
    1 分子について、指定された記述子を計算して dict で返す。
    mol が None の場合は、全部 None にして返す。
    """
    if mol is None:
        return {
            "MolWt": None,
            "NumAtoms": None,
            "NumHeavyAtoms": None,
            "NumHeteroAtoms": None,
            "NumBonds": None,
            "NumRotatableBonds": None,
            "NumHAcceptors": None,
            "NumHDonors": None,
            "BalabanJ": None,
            "BertzCT": None,
            "TPSA": None,
            "MolLogP": None,
        }

    # 個々の記述子を計算
    mol_wt = Descriptors.MolWt(mol)
    num_atoms = mol.GetNumAtoms()
    num_heavy = mol.GetNumHeavyAtoms()
    num_hetero = rdMolDescriptors.CalcNumHeteroatoms(mol)
    num_bonds = mol.GetNumBonds()
    num_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    num_h_acc = Lipinski.NumHAcceptors(mol)
    num_h_don = Lipinski.NumHDonors(mol)
    balaban_j = Descriptors.BalabanJ(mol)
    bertz_ct = Descriptors.BertzCT(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    mol_logp = Crippen.MolLogP(mol)

    return {
        "MolWt": mol_wt,
        "NumAtoms": num_atoms,
        "NumHeavyAtoms": num_heavy,
        "NumHeteroAtoms": num_hetero,
        "NumBonds": num_bonds,
        "NumRotatableBonds": num_rot,
        "NumHAcceptors": num_h_acc,
        "NumHDonors": num_h_don,
        "BalabanJ": balaban_j,
        "BertzCT": bertz_ct,
        "TPSA": tpsa,
        "MolLogP": mol_logp,
    }


# =========================
# 3. DataFrame 全体の組み立て
# =========================

DESCRIPTOR_COLUMNS = [
    "MolWt",
    "NumAtoms",
    "NumHeavyAtoms",
    "NumHeteroAtoms",
    "NumBonds",
    "NumRotatableBonds",
    "NumHAcceptors",
    "NumHDonors",
    "BalabanJ",
    "BertzCT",
    "TPSA",
    "MolLogP",
]


def build_descriptors_df(df_input: pd.DataFrame) -> pd.DataFrame:
    """
    ID/SMILES を持つ DataFrame を受け取り、
    RDKit 記述子列を追加した DataFrame を返す。
    """
    df_out = df_input.copy()

    # 記述子列を初期化
    for col in DESCRIPTOR_COLUMNS:
        if col not in df_out.columns:
            df_out[col] = None

    # 各行ごとに逐次処理
    for idx, row in df_out.iterrows():
        smiles = row["SMILES"]
        mol = smiles_to_mol(smiles)
        desc = calc_descriptors_for_mol(mol)

        for col in DESCRIPTOR_COLUMNS:
            df_out.at[idx, col] = desc.get(col)

    return df_out


def save_descriptors_excel(df: pd.DataFrame, out_path: Path) -> Path:
    """
    記述子付き DataFrame を out_path に Excel で保存して、
    そのパスを返す。
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(out_path, index=False)
    return out_path


# =========================
# 4. 外部公開関数（app.py から呼ぶ用）
# =========================


def run_descriptors_for_excel(progress_path: Path, result_dir: Path) -> Path:
    """
    app.py から呼び出すためのラッパ。

    Parameters
    ----------
    progress_path : Path
        progress.xlsx のパス（A列: ID, B列: SMILES を含む）
    result_dir : Path
        result フォルダ（例: Path("./result")）

    Returns
    -------
    Path
        result/progress_rdkit_descriptors.xlsx のパス
    """
    progress_path = progress_path.resolve()
    result_dir = result_dir.resolve()

    if not progress_path.is_file():
        raise FileNotFoundError(f"progress Excel が見つかりません: {progress_path}")

    # A,B列だけを読み込み（ID, SMILES）
    df_input = load_id_smiles_from_excel(progress_path)

    # RDKit 記述子を付与
    df_desc = build_descriptors_df(df_input)

    # result/progress_rdkit_descriptors.xlsx として保存
    out_path = make_output_path(result_dir)
    save_descriptors_excel(df_desc, out_path)

    return out_path
