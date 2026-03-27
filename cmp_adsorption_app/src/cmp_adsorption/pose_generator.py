from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase import Atoms


@dataclass
class PoseConfig:
    xy_step: float
    z_values: list[float]
    rotations_deg: list[float]
    max_poses: int


def _center_xy(atoms: Atoms) -> Atoms:
    out = atoms.copy()
    pos = out.get_positions()
    pos[:, :2] -= pos[:, :2].mean(axis=0)
    out.set_positions(pos)
    return out


def generate_initial_poses(slab: Atoms, adsorbate: Atoms, cfg: PoseConfig) -> list[Atoms]:
    slab_pos = slab.get_positions()
    min_xy = slab_pos[:, :2].min(axis=0)
    max_xy = slab_pos[:, :2].max(axis=0)
    xgrid = np.arange(min_xy[0], max_xy[0] + 1e-9, cfg.xy_step)
    ygrid = np.arange(min_xy[1], max_xy[1] + 1e-9, cfg.xy_step)
    z_top = slab_pos[:, 2].max()

    base = _center_xy(adsorbate)
    poses: list[Atoms] = []
    for x in xgrid:
        for y in ygrid:
            for z in cfg.z_values:
                for r in cfg.rotations_deg:
                    ads = base.copy()
                    ads.rotate(r, "z", center="COP")
                    apos = ads.get_positions()
                    apos += np.array([x, y, z_top + z])
                    ads.set_positions(apos)
                    combined = slab.copy() + ads
                    poses.append(combined)
                    if len(poses) >= cfg.max_poses:
                        return poses
    return poses
