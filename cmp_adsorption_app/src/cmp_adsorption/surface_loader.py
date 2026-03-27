from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import read


@dataclass
class SurfaceModel:
    sid: str
    slab: Atoms
    fixed_indices: list[int]


def _read_fixed_indices(surface_dir: Path) -> list[int]:
    candidates = [surface_dir / "fixed_index.txt", surface_dir / "fixed_indices.txt"]
    for p in candidates:
        if p.exists():
            values = []
            for line in p.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line:
                    continue
                values.append(int(line))
            return values
    return []


def load_surface(sid: str, surface_dir: Path) -> SurfaceModel:
    slab_path = surface_dir / "slab.xyz"
    slab = read(str(slab_path))
    fixed = _read_fixed_indices(surface_dir)
    if fixed:
        mask = np.zeros(len(slab), dtype=bool)
        mask[fixed] = True
        slab.set_constraint(FixAtoms(mask=mask))
    return SurfaceModel(sid=sid, slab=slab, fixed_indices=fixed)
