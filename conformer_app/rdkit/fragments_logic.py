#!/usr/bin/env python
from __future__ import annotations

"""
RDKit 官能基フラグメントカウント

指定 Excel ファイル (A列: ID, B列: SMILES) から
RDKit Chem.Fragments の fr_* フラグメント個数を計算し、
指定ディレクトリに ID/SMILES + fr_* 列を持つ Excel を出力する。

主なエントリポイント:
    run_fragments_for_excel(input_path: Path, output_dir: Path) -> Path

CLI としても使える:
    python rdkit_fragments_logic.py input.xlsx [--outdir DIR]
"""

import sys
import inspect
from pathlib import Path
from typing import List, Dict

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments


# ------------------------------
# 1. 入出力まわり
# ------------------------------


def load_id_smiles_from_excel(path: Path) -> pd.DataFrame:
    """
    Excel から A列(0), B列(1)だけを読み込んで DataFrame を返す。
    1行目はヘッダー、2行目以降がデータを想定。
    """
    df = pd.read_excel(path, usecols=[0, 1])
    # カラム名を固定しておく（元が何でも ID / SMILES にする）
    df.columns = ["ID", "SMILES"]
    return df


def make_output_filename(input_path: Path, output_dir: Path) -> Path:
    """
    入力ファイル名に _fragments を付けた Excel ファイル名を、
    指定された output_dir に作る。
        例: input_path = root/result/Progress.xlsx
            output_dir = root/result
            → root/result/Progress_fragments.xlsx
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / f"{input_path.stem}_rdkit_fragments.xlsx"


# ------------------------------
# 2. RDKit: フラグメント関数の準備
# ------------------------------


def get_fragment_functions():
    """
    rdkit.Chem.Fragments の中から fr_ で始まる関数を全部拾って
    [(名前, 関数), ...] のリストとして返す。
    """
    frag_funcs = [
        (name, func)
        for name, func in inspect.getmembers(Fragments, inspect.isfunction)
        if name.startswith("fr_")
    ]
    return frag_funcs


FRAG_FUNCS = get_fragment_functions()
FRAG_NAMES: List[str] = [name for name, _ in FRAG_FUNCS]


# ------------------------------
# 3. 1分子ぶんの計算
# ------------------------------


def smiles_to_mol(smiles: str):
    """SMILES 文字列から RDKit Mol を作る。失敗したら None を返す。"""
    if pd.isna(smiles):
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    return mol


def calc_fragments_for_mol(mol) -> Dict[str, int | None]:
    """
    1 分子について、全 fr_* の {名前: 個数} dict を返す。
    mol が None の場合は、全部 None にして返す。
    """
    if mol is None:
        return {name: None for name in FRAG_NAMES}

    return {name: func(mol) for name, func in FRAG_FUNCS}


# ------------------------------
# 4. 行ごとに逐次処理して DF を作る部分
# ------------------------------


def build_fragments_dataframe(input_path: Path) -> pd.DataFrame:
    """
    入力 Excel から ID/SMILES を読み込み、
    fr_* 列を追加した DataFrame を返す（まだファイルには書き出さない）。
    """
    # A,B列だけ読み込む
    df_input = load_id_smiles_from_excel(input_path)

    # 出力用 DataFrame: とりあえず ID/SMILES だけコピー
    df_out = df_input.copy()

    # フラグメント用の列を追加（初期値は None）
    for name in FRAG_NAMES:
        df_out[name] = None

    # 各行ごとに逐次処理
    for idx, row in df_out.iterrows():
        smiles = row["SMILES"]
        mol = smiles_to_mol(smiles)
        frags = calc_fragments_for_mol(mol)

        # 結果をその行に書き込む
        for name in FRAG_NAMES:
            df_out.at[idx, name] = frags[name]

    return df_out


def run_fragments_for_excel(input_path: Path, output_dir: Path) -> Path:
    """
    GUI から呼び出す想定のメイン関数。

    Parameters
    ----------
    input_path : Path
        フラグメント解析の元になる Excel ファイル
        （例: root/result/Progress.xlsx）
    output_dir : Path
        出力先ディレクトリ
        （例: root/result）

    Returns
    -------
    out_path : Path
        出力された Excel のパス
        （例: root/result/Progress_fragments.xlsx）
    """
    input_path = input_path.resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"入力ファイルが見つかりません: {input_path}")

    df_out = build_fragments_dataframe(input_path)
    out_path = make_output_filename(input_path, output_dir)
    df_out.to_excel(out_path, index=False)
    print(f"[fragments] 出力ファイル: {out_path}")
    return out_path


# ------------------------------
# 5. CLI (オプション)
# ------------------------------


def main():
    """
    CLI 用エントリポイント。
    デフォルトでは input.xlsx と同じディレクトリに *_fragments.xlsx を出力。
    """
    if len(sys.argv) < 2:
        print("使い方: python rdkit_fragments_logic.py input.xlsx [outdir]")
        sys.exit(1)

    input_path = Path(sys.argv[1])
    if len(sys.argv) >= 3:
        output_dir = Path(sys.argv[2])
    else:
        # デフォルトは input.xlsx と同じディレクトリ
        output_dir = input_path.parent

    try:
        run_fragments_for_excel(input_path, output_dir)
    except Exception as e:
        print(f"エラー: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
