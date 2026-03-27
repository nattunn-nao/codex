from __future__ import annotations

from ase import Atoms


def min_interatomic_distance(structure: Atoms, slab_size: int) -> float:
    pos = structure.get_positions()
    m = 1e9
    for i in range(slab_size):
        for j in range(slab_size, len(structure)):
            d = ((pos[i] - pos[j]) ** 2).sum() ** 0.5
            if d < m:
                m = d
    return float(m)


def geometry_pass(structure: Atoms, slab_size: int, min_distance: float) -> bool:
    return min_interatomic_distance(structure, slab_size) >= min_distance
