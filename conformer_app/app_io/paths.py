# conformer_app/app_io/paths.py
from __future__ import annotations

from pathlib import Path
from typing import Union

PathLike = Union[str, Path]


def _ensure_dir(path: PathLike) -> Path:
    """
    パスに対応するディレクトリを作成（親もまとめて）し、
    Path オブジェクトを返すユーティリティ。
    """
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


# =========================
# ID ごとのルート
# =========================

def id_root_dir(base_dir: PathLike, mol_id: str) -> Path:
    """
    output ベースディレクトリと ID から
    <base_dir>/<ID>/ を作成して返す。
    例: ./output/111-200
    """
    return _ensure_dir(Path(base_dir) / str(mol_id))


# =========================
# サブディレクトリ群（ID 配下）
# =========================

def logs_dir(id_dir: PathLike) -> Path:
    """<ID>/logs/"""
    return _ensure_dir(Path(id_dir) / "logs")


def mm_dir(id_dir: PathLike) -> Path:
    """<ID>/MM/"""
    return _ensure_dir(Path(id_dir) / "MM")


def post_mm_dir(id_dir: PathLike) -> Path:
    """<ID>/post_MM/  … 重複除去・エネルギーソート後の構造置き場"""
    return _ensure_dir(Path(id_dir) / "post_MM")


def uma_dir(id_dir: PathLike) -> Path:
    """<ID>/UMA/"""
    return _ensure_dir(Path(id_dir) / "UMA")


def qm_dir(id_dir: PathLike) -> Path:
    """<ID>/QM/"""
    return _ensure_dir(Path(id_dir) / "QM")


# =========================
# result 用
# =========================

def result_root(base_dir: PathLike) -> Path:
    """
    result のベースディレクトリを作成して返す。
    例: ./result
    """
    return _ensure_dir(Path(base_dir))


def result_root_dir(base_dir: PathLike) -> Path:
    """
    互換用エイリアス。
    以前のコードで使っていた result_root_dir を
    result_root にフォワードする。
    """
    return result_root(base_dir)


def result_dir(base_dir: PathLike, mol_id: str) -> Path:
    """
    result ベースディレクトリと ID から
    <base_dir>/<ID>/ を作成して返す。
    例: ./result/111-200
    """
    root = result_root(base_dir)
    return _ensure_dir(root / str(mol_id))
