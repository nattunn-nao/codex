# conformer_app/qm_layer/input_builder.py
from __future__ import annotations

from pathlib import Path
from typing import List, Sequence

from ase.io import read


# ============================================================
# 基本ユーティリティ
# ============================================================


def _find_project_root(start: Path) -> Path:
    """
    app.py がある階層、または templates ディレクトリがある階層を
    「プロジェクトルート」とみなして探索する。
    例:
        <root>/app.py
        <root>/templates/gms.inp
        <root>/conformer_app/qm_layer/input_builder.py
    """
    cur = start
    for _ in range(6):  # 6 階層くらいまで上にたどる
        if (cur / "app.py").is_file() or (cur / "templates").is_dir():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return start


def _resolve_template_path(tpl_hint: str) -> str:
    """
    gms.inp テンプレートのパスを解決する。

    優先順位:
      1) ユーザー指定パス (tpl_hint) またはその .inp 付き
      2) プロジェクトルートの templates/gms.inp   ← ★ app.py 基準 ./templates/gms.inp
      3) プロジェクトルートの template/gms.inp    （古い構成との互換用）
      4) カレントディレクトリの templates/gms.inp / template/gms.inp
    """
    candidates: List[Path] = []

    # 1) UI 等で指定されたヒント（相対パスでも絶対パスでもOK）
    if tpl_hint:
        p = Path(tpl_hint).expanduser()
        candidates.append(p)
        if p.suffix == "":
            candidates.append(p.with_suffix(".inp"))

    # 2) 3) プロジェクトルートを推定してそこから探す
    here = Path(__file__).resolve()
    # 例: <root>/conformer_app/qm_layer/input_builder.py
    #   → here.parent = qm_layer
    project_root = _find_project_root(here.parent)

    # app.py と同じ階層にある ./templates/gms.inp を最優先で見る
    candidates.append(project_root / "templates" / "gms.inp")  # ★あなたの配置
    candidates.append(project_root / "template" / "gms.inp")  # 互換用

    # 4) 実行時カレントディレクトリ基準でも探す（streamlit run の場所がズレている場合用）
    cwd = Path.cwd()
    candidates.append(cwd / "templates" / "gms.inp")
    candidates.append(cwd / "template" / "gms.inp")

    for c in candidates:
        if c.is_file():
            return str(c)

    msg_lines = ["gms.inp テンプレートが見つかりませんでした。探索した候補:"]
    msg_lines += [f"  - {c}" for c in candidates]
    raise FileNotFoundError("\n".join(msg_lines))


def _read_xyz_first(path: str):
    """
    ASE の read() を使って XYZ を読み、最初の構造だけを返す。
    （マルチフレーム XYZ だった場合の安全策）
    """
    obj = read(path, format="xyz", index=":")
    if isinstance(obj, list):
        return obj[0]
    return obj


def xyz_to_gamess_unique_block(xyz_path: str) -> str:
    """
    XYZ ファイルから GAMESS の COORD=UNIQUE 用ブロックを生成する。

    形式:
        Sym   Z    x   y   z
    のような行を連結した文字列を返す。
    """
    atoms = _read_xyz_first(xyz_path)
    zs = atoms.get_atomic_numbers()
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()

    lines: List[str] = []
    for i, (sym, pos) in enumerate(zip(symbols, positions)):
        Z = float(zs[i])
        x, y, z = pos
        lines.append(f"{sym:<2s}  {Z:.1f}    {x: .8f}  {y: .8f}  {z: .8f}")
    return "\n".join(lines)


def fill_gamess_template(
    template_text: str,
    comment: str,
    pointgroup: str,
    gen_xyz_block: str,
    icharge: int,
    mult: int,
) -> str:
    """
    gms.inp のテンプレートテキスト中のプレースホルダを実際の値で置換する。

    想定プレースホルダ:
      - %COMMENT_SMILES_N%
      - %POINTGROUP%
      - %GEN_XYZ%
      - %ICHARG%
      - %MULT%
    """
    out = (
        template_text
        .replace("%COMMENT_SMILES_N%", comment)
        .replace("%POINTGROUP%", pointgroup)
        .replace("%GEN_XYZ%", gen_xyz_block)
        .replace("%ICHARG%", str(int(icharge)))
        .replace("%MULT%", str(int(mult)))
    )
    return out


# ============================================================
# .inp 生成メイン関数
# ============================================================


def make_inp_for_slots(
    id_dir: str | Path,
    gms_template_path: str,
    smiles_comment: str = "",
    pointgroup: str = "C1",
    icharge: int = 0,
    mult: int = 1,
    slots: Sequence[str] = ("gms1", "gms2", "gms3"),
    geom_filename: str = "geom.xyz",
) -> List[str]:
    """
    1つの ID に対して、QM/gms1..gmsN フォルダ内の geom.xyz から
    GAMESS 用の .inp を生成する。

    Parameters
    ----------
    id_dir : str or Path
        1つの ID に対応する output フォルダ
        例: ./output/111-1
    gms_template_path : str
        gms.inp テンプレートへのヒントパス。
        - "" や "./template/gms.inp" などが渡される想定。
        - 実際には app.py と同じ階層の ./templates/gms.inp を最優先で探索する。
    smiles_comment : str
        テンプレート内 %COMMENT_SMILES_N% に埋め込む文字列。
        空文字の場合は id_dir の basename を使うとよい。
    pointgroup : str
        %POINTGROUP% に入れるポイントグループ文字列（例: "C1"）。
    icharge : int
        GAMESS の ICHARG。
    mult : int
        GAMESS の MULT。
    slots : Sequence[str]
        処理対象となる gms サブディレクトリ名のリスト。
        例: ("gms1", "gms2", "gms3")
    geom_filename : str
        QM/gmsX/ 配下で座標として読むファイル名（デフォルト "geom.xyz"）。

    Returns
    -------
    out_paths : list of str
        生成された .inp ファイルのフルパスリスト。
        1つも生成できなかった場合は FileNotFoundError を送出する。
    """
    id_dir = Path(id_dir)
    qm_dir = id_dir / "QM"

    # テンプレートパスの解決（UI指定 or プロジェクトルート ./templates/gms.inp）
    tpl_path = Path(_resolve_template_path(gms_template_path))

    if not tpl_path.is_file():
        raise FileNotFoundError(f"gms.inp テンプレートが見つかりませんでした: {tpl_path}")

    template_text = tpl_path.read_text(encoding="utf-8")

    # コメントが空なら ID 名をデフォルトで入れる
    if not smiles_comment:
        smiles_comment = id_dir.name

    out_paths: List[str] = []

    for slot in slots:
        sdir = qm_dir / slot
        xyz = sdir / geom_filename
        if not xyz.is_file():
            # このスロットには geom.xyz が無ければスキップ
            continue

        gen_xyz_block = xyz_to_gamess_unique_block(str(xyz))
        inp_text = fill_gamess_template(
            template_text=template_text,
            comment=smiles_comment,
            pointgroup=pointgroup,
            gen_xyz_block=gen_xyz_block,
            icharge=icharge,
            mult=mult,
        )

        sdir.mkdir(parents=True, exist_ok=True)
        outp = sdir / f"{slot}.inp"
        outp.write_text(inp_text, encoding="utf-8")
        out_paths.append(str(outp))

    if not out_paths:
        raise FileNotFoundError(
            f"{qm_dir} 配下に {geom_filename} が見つからず、.inp を生成できませんでした。"
        )

    return out_paths
