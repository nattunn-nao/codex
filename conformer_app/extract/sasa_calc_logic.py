from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
from skimage.measure import marching_cubes

# SASA 用プローブ半径 (Å)
PROBE_RADIUS: float = 1.4

# 距離場グリッドの分解能 (Å)
GRID_SPACING: float = 0.25


# ============================================================
# グリッド & 距離場
# ============================================================

def make_grid(coords: np.ndarray,
              radii: np.ndarray,
              spacing: float = GRID_SPACING
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                         np.ndarray, np.ndarray, np.ndarray]:
    """
    全原子を覆う直方体グリッドを作成。
    戻り値: xs, ys, zs, X, Y, Z
    """
    coords = np.asarray(coords, float)
    max_r = float(radii.max())
    margin = max_r + 0.5

    mins = coords.min(axis=0) - margin
    maxs = coords.max(axis=0) + margin

    xs = np.arange(mins[0], maxs[0], spacing)
    ys = np.arange(mins[1], maxs[1], spacing)
    zs = np.arange(mins[2], maxs[2], spacing)

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    return xs, ys, zs, X, Y, Z


def distance_field_on_grid(X: np.ndarray,
                           Y: np.ndarray,
                           Z: np.ndarray,
                           coords: np.ndarray,
                           radii: np.ndarray) -> np.ndarray:
    """
    球の union に対する signed distance field を計算。
      field < 0 : 球の内側
      field = 0 : 球の表面
      field > 0 : 球の外側
    """
    if len(coords) == 0:
        return np.ones_like(X, dtype=float)

    coords = np.asarray(coords, float)
    radii = np.asarray(radii, float)

    grid_points = np.stack([X, Y, Z], axis=-1)  # (nx,ny,nz,3)
    diff = grid_points[..., None, :] - coords[None, None, None, :, :]  # (...,N,3)
    dist = np.linalg.norm(diff, axis=-1)  # (nx,ny,nz,N)
    d = dist - radii                       # signed distance
    field = np.min(d, axis=-1)            # union: min
    return field


# ============================================================
# メッシュ生成
# ============================================================

def field_to_vertices_faces(field: np.ndarray,
                            xs: np.ndarray,
                            ys: np.ndarray,
                            zs: np.ndarray,
                            level: float = 0.0,
                            step_size: int = 1
                            ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    signed distance field から level=0 の等値面メッシュを生成。
    戻り値:
      verts : (n_verts, 3)
      faces : (n_faces, 3) int
    """
    vmin, vmax = float(field.min()), float(field.max())
    if not (vmin <= level <= vmax):
        return None, None

    dx = xs[1] - xs[0] if len(xs) > 1 else 1.0
    dy = ys[1] - ys[0] if len(ys) > 1 else 1.0
    dz = zs[1] - zs[0] if len(zs) > 1 else 1.0

    verts, faces, _, _ = marching_cubes(
        field,
        level=level,
        spacing=(dx, dy, dz),
        step_size=step_size,
    )

    verts[:, 0] += xs[0]
    verts[:, 1] += ys[0]
    verts[:, 2] += zs[0]

    return verts, faces.astype(int)


# ============================================================
# メッシュから面積・体積など
# ============================================================

def mesh_area(verts: Optional[np.ndarray],
              faces: Optional[np.ndarray]) -> float:
    """三角形メッシュの総表面積 (Å^2) を返す"""
    if verts is None or faces is None:
        return np.nan
    tri = verts[faces]                  # (n_faces, 3, 3)
    v1 = tri[:, 1] - tri[:, 0]
    v2 = tri[:, 2] - tri[:, 0]
    area = 0.5 * np.linalg.norm(np.cross(v1, v2), axis=1).sum()
    return float(area)


def mesh_area_polar(verts: Optional[np.ndarray],
                    faces: Optional[np.ndarray],
                    vertex_polar: Optional[np.ndarray]) -> float:
    """
    極性部 SASA 面積 (Å^2) を返す。
    各三角形について、頂点いずれかが極性ならその三角形を極性領域とみなす。
    """
    if verts is None or faces is None or vertex_polar is None:
        return np.nan

    tri = verts[faces]  # (n_faces, 3, 3)
    v1 = tri[:, 1] - tri[:, 0]
    v2 = tri[:, 2] - tri[:, 0]
    area_each = 0.5 * np.linalg.norm(np.cross(v1, v2), axis=1)  # (n_faces,)

    face_polar = vertex_polar[faces].any(axis=1)  # (n_faces,)
    polar_area = area_each[face_polar].sum()
    return float(polar_area)


def volume_from_field(field: Optional[np.ndarray],
                      xs: np.ndarray,
                      ys: np.ndarray,
                      zs: np.ndarray) -> float:
    """field < 0 の体積 (Å^3) を粗くボクセルから見積もる"""
    if field is None:
        return np.nan
    dx = xs[1] - xs[0] if len(xs) > 1 else 1.0
    dy = ys[1] - ys[0] if len(ys) > 1 else 1.0
    dz = zs[1] - zs[0] if len(zs) > 1 else 1.0
    voxel_volume = dx * dy * dz
    inside = (field < 0.0)
    return float(inside.sum() * voxel_volume)


# ============================================================
# 頂点の極性判定（SASA 計算に使用）
# ============================================================

def classify_vertices_polar(verts: Optional[np.ndarray],
                            atom_coords: np.ndarray,
                            polar_mask: np.ndarray) -> Optional[np.ndarray]:
    """
    各頂点について、最も近い原子が極性かどうかを判定。
    戻り値: vertex_polar (bool 配列, len = n_verts)
    """
    if verts is None:
        return None

    atom_coords = np.asarray(atom_coords, float)
    diff = verts[:, None, :] - atom_coords[None, :, :]  # (n_verts, n_atoms, 3)
    dist2 = np.sum(diff**2, axis=-1)                    # (n_verts, n_atoms)
    nearest_atom_idx = np.argmin(dist2, axis=1)         # (n_verts,)
    return polar_mask[nearest_atom_idx]


# ============================================================
# SASA + vdW の計算をまとめた高レベル関数
# ============================================================

def compute_sasa_vdw_properties(coords: np.ndarray,
                                polar_mask: np.ndarray,
                                xs: np.ndarray,
                                ys: np.ndarray,
                                zs: np.ndarray,
                                field_sasa: np.ndarray,
                                field_vdw: np.ndarray) -> Dict[str, object]:
    """
    既に計算済みの field_sasa / field_vdw とグリッドから、
    - SASA 表面メッシュ
    - vdW 表面メッシュ
    - SASA 面積、極性 SASA 面積
    - vdW 面積、vdW 体積
    などをまとめて返す。
    """
    verts_sasa, faces_sasa = field_to_vertices_faces(
        field_sasa, xs, ys, zs, level=0.0, step_size=1
    )
    verts_vdw, faces_vdw = field_to_vertices_faces(
        field_vdw, xs, ys, zs, level=0.0, step_size=1
    )

    vertex_polar_sasa = classify_vertices_polar(verts_sasa, coords, polar_mask)

    sasa_area = mesh_area(verts_sasa, faces_sasa)
    sasa_polar_area = mesh_area_polar(verts_sasa, faces_sasa, vertex_polar_sasa)
    vdw_area = mesh_area(verts_vdw, faces_vdw)
    vdw_volume = volume_from_field(field_vdw, xs, ys, zs)

    return {
        "verts_sasa": verts_sasa,
        "faces_sasa": faces_sasa,
        "verts_vdw": verts_vdw,
        "faces_vdw": faces_vdw,
        "vertex_polar_sasa": vertex_polar_sasa,
        "sasa_area": sasa_area,
        "sasa_polar_area": sasa_polar_area,
        "vdw_area": vdw_area,
        "vdw_volume": vdw_volume,
    }
