# conformer_app/app_io/excel_io.py
from __future__ import annotations

from pathlib import Path
from typing import Any, List, Sequence, Tuple, Union

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

# RDKit の警告を消しておく（ノイズ防止）
RDLogger.DisableLog("rdApp.*")


# ======================================================================
# 内部ユーティリティ
# ======================================================================

def _compute_formal_charge(smiles: str) -> int:
    """
    SMILES から形式電荷（formal charge）を計算する。

    - RDKit で SMILES をパースし、全原子の formal charge の総和を返す
    - パースに失敗した場合や何か例外が起きた場合は 0 を返す
    - 注意: SMILES に電荷が明示されていない場合（例: 'CN'）は 0 になる
      （pH によるプロトン化状態までは自動推定しない）
    """
    smiles = (smiles or "").strip()
    if not smiles:
        return 0

    mol = None
    try:
        # 通常パース
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # sanitize=False で再トライしてから Sanitize してみる
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                Chem.SanitizeMol(mol)
    except Exception:
        mol = None

    if mol is None:
        # パースに失敗したら中性扱い
        return 0

    try:
        # 分子全体の formal charge（原子の GetFormalCharge の総和）
        # RDKit の Mol.GetFormalCharge はこの総和を返す
        chg = int(mol.GetFormalCharge())
        return chg
    except Exception:
        # 念のためのフォールバック
        try:
            total = 0
            for atom in mol.GetAtoms():
                total += int(atom.GetFormalCharge())
            return total
        except Exception:
            return 0


def _parse_requirement_df(df: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """
    requirement 用の DataFrame から ID / SMILES リストを取り出す。
    A列=ID, B列=SMILES を想定。ヘッダー名は問わない。
    1行目が "SMILES" などのヘッダー行でも自動でスキップする。
    """
    if df.shape[1] < 2:
        raise ValueError("Excel の列数が足りません（少なくとも A列=ID, B列=SMILES が必要）。")

    ids_raw = df.iloc[:, 0].astype(str).str.strip().tolist()
    smis_raw = df.iloc[:, 1].astype(str).str.strip().tolist()

    ids: List[str] = []
    smis: List[str] = []

    for id_raw, smi_raw in zip(ids_raw, smis_raw):
        if not smi_raw or smi_raw.lower() == "smiles":
            # SMILES 列のヘッダ行や空行はスキップ
            continue
        if not id_raw:
            id_raw = smi_raw
        ids.append(id_raw)
        smis.append(smi_raw)

    return ids, smis


# ======================================================================
# requirement Excel 読み込み
# ======================================================================

def parse_requirement_excel(uploaded_file) -> Tuple[List[str], List[str]]:
    """
    requirement Excel (.xlsx) から ID と SMILES のリストを取り出す。

    想定フォーマット:
      - A列: ID
      - B列: SMILES
      - C列以降: 無視

    ヘッダー行の有無は問わず、"SMILES" という文字列を含む行はスキップする。
    """
    if uploaded_file is None:
        return [], []

    try:
        df = pd.read_excel(uploaded_file)
    except Exception as e:
        raise RuntimeError(f"Excel の読み込みに失敗しました: {e}") from e

    ids, smis = _parse_requirement_df(df)
    return ids, smis


def load_requirement_excel(path: str | Path) -> Tuple[List[str], List[str]]:
    """
    requirement Excel (.xlsx) を読み込み、(ids, smis) を返す。

    - A列: ID
    - B列: SMILES
    - C列以降は無視
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"requirement Excel が見つかりません: {path}")

    try:
        df = pd.read_excel(path)
    except Exception as e:
        raise RuntimeError(f"requirement Excel の読み込みに失敗しました: {path} : {e}") from e

    ids, smis = _parse_requirement_df(df)
    if not ids:
        raise RuntimeError(f"requirement Excel に有効な ID/SMILES 行が見つかりませんでした: {path}")

    return ids, smis


# ======================================================================
# 初期 Progress.xlsx 作成
# ======================================================================

def create_initial_progress_excel(
    ids: Sequence[str],
    smis: Sequence[str],
    out_path: str | Path,
) -> str:
    """
    requirement Excel から読み込んだ ID / SMILES から、
    計算開始前に progress.xlsx の「雛形」を作成する。

    カラム構成:
      - ID
      - SMILES
      - Formal_Charge : SMILES から計算した形式電荷（整数, 失敗時は 0）
      - MM_status
      - UMA_status
      - QM_status
      - error_message
    """
    records: List[dict[str, Any]] = []

    for id_raw, smi_raw in zip(ids, smis):
        id_str = str(id_raw).strip()
        smi_str = str(smi_raw).strip()
        if not smi_str:
            continue

        chg = _compute_formal_charge(smi_str)

        records.append(
            {
                "ID": id_str,
                "SMILES": smi_str,
                "Formal_Charge": chg,
                "MM_status": "WAIT",
                "UMA_status": "WAIT",
                "QM_status": "WAIT",
                "error_message": "",
            }
        )

    df = pd.DataFrame(records, columns=[
        "ID",
        "SMILES",
        "Formal_Charge",
        "MM_status",
        "UMA_status",
        "QM_status",
        "error_message",
    ])

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(out_path, index=False, sheet_name="progress")

    return str(out_path)


# ======================================================================
# Progress.xlsx 読み込み（再開機能などで利用）
# ======================================================================

def load_progress_excel(path: str | Path) -> pd.DataFrame:
    """
    既存の progress.xlsx を読み込んで DataFrame を返す。
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"progress.xlsx が見つかりません: {path}")
    try:
        return pd.read_excel(path)
    except Exception as e:
        raise RuntimeError(f"progress.xlsx の読み込みに失敗しました: {path} : {e}") from e


# ======================================================================
# 実行中 Progress.xlsx 更新（部分的上書き）
# ======================================================================

def write_progress_excel(
    rows: Sequence[Any],
    out_path: str | Path,
) -> str:
    """
    JobProgressRow 相当のオブジェクト列から progress.xlsx を更新する。

    仕様:
      - 既存の progress.xlsx を読み込み
      - rows に含まれる ID と一致する行だけ、
        MM_status / UMA_status / QM_status / error_message / Formal_Charge を更新
      - それ以外の ID の行はそのまま残す（= 未処理 ID 行は消えない）

    rows の各要素は少なくとも以下の属性を持つことを想定:
      - id
      - smiles
      - mm_status
      - uma_status
      - qm_status
      - error_message

    mm_status / uma_status / qm_status は Enum(StageState) でも文字列でもよい。
    """
    out_path = Path(out_path)

    # progress.xlsx が存在しない場合は、rows だけで新規作成（最後の手段）
    if not out_path.is_file():
        records: List[dict[str, Any]] = []

        for r in rows:
            mol_id = getattr(r, "id", "")
            smi = getattr(r, "smiles", "")
            chg = _compute_formal_charge(str(smi))

            def _to_state_str(x: Any) -> str:
                if hasattr(x, "value"):
                    try:
                        return str(x.value)
                    except Exception:
                        return str(x)
                return str(x)

            mm_status = _to_state_str(getattr(r, "mm_status", ""))
            uma_status = _to_state_str(getattr(r, "uma_status", ""))
            qm_status = _to_state_str(getattr(r, "qm_status", ""))
            err = getattr(r, "error_message", "")

            records.append(
                {
                    "ID": str(mol_id),
                    "SMILES": str(smi),
                    "Formal_Charge": chg,
                    "MM_status": mm_status,
                    "UMA_status": uma_status,
                    "QM_status": qm_status,
                    "error_message": str(err) if err is not None else "",
                }
            )

        df_new = pd.DataFrame(records, columns=[
            "ID",
            "SMILES",
            "Formal_Charge",
            "MM_status",
            "UMA_status",
            "QM_status",
            "error_message",
        ])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df_new.to_excel(out_path, index=False, sheet_name="progress")
        return str(out_path)

    # ここからが通常ルート：既存 progress.xlsx を部分更新
    df = load_progress_excel(out_path)

    # ID 列が無ければ fallback として rows だけで再作成
    if "ID" not in df.columns:
        records: List[dict[str, Any]] = []
        for r in rows:
            mol_id = getattr(r, "id", "")
            smi = getattr(r, "smiles", "")
            chg = _compute_formal_charge(str(smi))

            def _to_state_str(x: Any) -> str:
                if hasattr(x, "value"):
                    try:
                        return str(x.value)
                    except Exception:
                        return str(x)
                return str(x)

            mm_status = _to_state_str(getattr(r, "mm_status", ""))
            uma_status = _to_state_str(getattr(r, "uma_status", ""))
            qm_status = _to_state_str(getattr(r, "qm_status", ""))
            err = getattr(r, "error_message", "")

            records.append(
                {
                    "ID": str(mol_id),
                    "SMILES": str(smi),
                    "Formal_Charge": chg,
                    "MM_status": mm_status,
                    "UMA_status": uma_status,
                    "QM_status": qm_status,
                    "error_message": str(err) if err is not None else "",
                }
            )
        df_new = pd.DataFrame(records, columns=[
            "ID",
            "SMILES",
            "Formal_Charge",
            "MM_status",
            "UMA_status",
            "QM_status",
            "error_message",
        ])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df_new.to_excel(out_path, index=False, sheet_name="progress")
        return str(out_path)

    # 各行ごとに ID をキーとして既存 DataFrame を更新
    for r in rows:
        mol_id = str(getattr(r, "id", "")).strip()
        smi = str(getattr(r, "smiles", "")).strip()
        if not mol_id:
            continue

        # 対象行を検索
        mask = df["ID"].astype(str).str.strip() == mol_id
        if not mask.any():
            # 既存 progress に無い ID の場合は append
            chg = _compute_formal_charge(smi)

            def _to_state_str(x: Any) -> str:
                if hasattr(x, "value"):
                    try:
                        return str(x.value)
                    except Exception:
                        return str(x)
                return str(x)

            mm_status = _to_state_str(getattr(r, "mm_status", ""))
            uma_status = _to_state_str(getattr(r, "uma_status", ""))
            qm_status = _to_state_str(getattr(r, "qm_status", ""))
            err = getattr(r, "error_message", "")

            new_row = {
                "ID": mol_id,
                "SMILES": smi,
                "Formal_Charge": chg,
                "MM_status": mm_status,
                "UMA_status": uma_status,
                "QM_status": qm_status,
                "error_message": str(err) if err is not None else "",
            }
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
            continue

        # 既存行を更新（複数あれば全部）
        chg = _compute_formal_charge(smi)

        def _to_state_str(x: Any) -> str:
            if hasattr(x, "value"):
                try:
                    return str(x.value)
                except Exception:
                    return str(x)
            return str(x)

        mm_status = _to_state_str(getattr(r, "mm_status", ""))
        uma_status = _to_state_str(getattr(r, "uma_status", ""))
        qm_status = _to_state_str(getattr(r, "qm_status", ""))
        err = getattr(r, "error_message", "")

        if "SMILES" in df.columns:
            df.loc[mask, "SMILES"] = smi
        if "Formal_Charge" in df.columns:
            df.loc[mask, "Formal_Charge"] = chg
        if "MM_status" in df.columns:
            df.loc[mask, "MM_status"] = mm_status
        if "UMA_status" in df.columns:
            df.loc[mask, "UMA_status"] = uma_status
        if "QM_status" in df.columns:
            df.loc[mask, "QM_status"] = qm_status
        if "error_message" in df.columns:
            df.loc[mask, "error_message"] = str(err) if err is not None else ""

    # 更新後の DataFrame を書き戻し
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(out_path, index=False, sheet_name="progress")

    return str(out_path)


# ======================================================================
# 再開用フラグ適用（apply_resume_flags_to_jobs）
# ======================================================================

def apply_resume_flags_to_jobs(
    jobs: Sequence[Any],
    resume_source: Union[str, Path, pd.DataFrame, None],
) -> List[Any]:
    """
    Progress / resume 情報をもとに「すでに完了した ID をスキップする」ためのヘルパ。

    jobs:
        - JobProgressRow のようなオブジェクト列
          （少なくとも .id を持つ）
        - または (id, smiles, ...) のようなタプル列
          （先頭要素を ID とみなす）

    resume_source:
        - progress.xlsx のパス (str / Path)
        - 既に読み込まれた DataFrame
        - None の場合は jobs をそのまま返す

    動作:
        - resume_source から "ID" と "QM_status" を読み取り、
          QM_status が "DONE"（大文字小文字無視）の ID はスキップ対象とする。
        - 上記以外（ERROR / ABORTED / 空欄など）は「未完了」と見なして実行対象とする。
        - 返り値は「スキップされなかった jobs のリスト」。
    """
    if resume_source is None:
        return list(jobs)

    # resume DataFrame を取得
    if isinstance(resume_source, pd.DataFrame):
        df = resume_source.copy()
    else:
        # パスとして扱う
        path = Path(resume_source)
        if not path.is_file():
            # 指定されているがファイルが無ければ、そのまま jobs を返す
            return list(jobs)
        df = load_progress_excel(path)

    cols = set(df.columns)
    if "ID" not in cols or "QM_status" not in cols:
        # 必要な列が無い場合は何もしない
        return list(jobs)

    # 完了済み ID の集合（QM_status == DONE）
    completed_ids = set(
        str(id_).strip()
        for id_, st in zip(df["ID"], df["QM_status"])
        if isinstance(st, str) and st.strip().upper() == "DONE"
    )

    if not completed_ids:
        return list(jobs)

    filtered: List[Any] = []

    for job in jobs:
        # job から ID を取り出す
        if hasattr(job, "id"):
            jid = str(getattr(job, "id"))
        elif isinstance(job, (tuple, list)) and len(job) > 0:
            jid = str(job[0])
        else:
            # ID がわからないものは念のため実行対象に残す
            filtered.append(job)
            continue

        if jid in completed_ids:
            # この ID は既に QM DONE と判断し、スキップする
            continue

        filtered.append(job)

    return filtered
