# app_io/xyz_io.py
from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
from ase import Atoms
from ase.io import read, write


PathLike = Union[str, Path]


def read_xyz_to_atoms(path: PathLike) -> Atoms:
    """
    XYZ ファイルを ASE Atoms として読み込む。
    """
    return read(str(path))


def _move_to_com(atoms: Atoms) -> None:
    """
    原子座標を重心が原点にくるように平行移動する（in-place）。
    """
    com = atoms.get_center_of_mass()
    pos = atoms.get_positions()
    atoms.set_positions(pos - com)


def write_atoms_to_xyz(
    atoms: Atoms,
    path: PathLike,
    comment: str = "",
    move_to_com: bool = False,
) -> None:
    """
    ASE Atoms を XYZ として書き出す。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    atoms_w = atoms.copy()
    if move_to_com:
        _move_to_com(atoms_w)

    write(str(path), atoms_w, format="xyz", comment=comment)
