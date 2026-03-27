from __future__ import annotations

from pathlib import Path

from ase import Atoms
from ase.io import read


def list_adsorbate_xyz(input_dir: Path) -> list[Path]:
    if not input_dir.exists():
        return []
    return sorted(input_dir.glob("*.xyz"))


def load_adsorbate(path: Path) -> Atoms:
    return read(str(path))
