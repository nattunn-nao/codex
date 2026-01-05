from __future__ import annotations

from typing import Optional, Tuple, List

import numpy as np
import plotly.graph_objects as go


# ワイヤーフレームの線の太さ
WIRE_LINE_WIDTH: int = 1


# ============================================================
# 色
# ============================================================

def element_color(sym: str) -> str:
    table = {
        "H": "white", "C": "gray", "N": "blue", "O": "red",
        "F": "green", "P": "orange", "S": "yellow",
        "Cl": "green", "Br": "brown", "I": "purple",
    }
    return table.get(sym, "lightgrey")


# ============================================================
# メッシュ → ワイヤーフレーム線分（極性/非極性で色分け）
# ============================================================

def faces_to_wireframe_lines_colored(verts: Optional[np.ndarray],
                                     faces: Optional[np.ndarray],
                                     vertex_polar: Optional[np.ndarray]
                                     ) -> Tuple[
                                         Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray],
                                         Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]
                                     ]:
    """
    三角形メッシュ (verts, faces) から
    極性/非極性エッジに分けたワイヤーフレーム座標列を作る。
    """
    if verts is None or faces is None or vertex_polar is None:
        return None, None, None, None, None, None

    edges = set()
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            u, v = sorted((int(u), int(v)))
            edges.add((u, v))

    x_pol, y_pol, z_pol = [], [], []
    x_non, y_non, z_non = [], [], []

    for u, v in edges:
        is_polar = bool(vertex_polar[u] or vertex_polar[v])
        xs = [verts[u, 0], verts[v, 0], None]
        ys = [verts[u, 1], verts[v, 1], None]
        zs = [verts[u, 2], verts[v, 2], None]

        if is_polar:
            x_pol.extend(xs)
            y_pol.extend(ys)
            z_pol.extend(zs)
        else:
            x_non.extend(xs)
            y_non.extend(ys)
            z_non.extend(zs)

    return (
        np.array(x_pol), np.array(y_pol), np.array(z_pol),
        np.array(x_non), np.array(y_non), np.array(z_non),
    )


# ============================================================
# プロット本体
# ============================================================

def plot_wireframe_colored(comment: str,
                           symbols: List[str],
                           coords: np.ndarray,
                           verts_sasa: Optional[np.ndarray],
                           faces_sasa: Optional[np.ndarray],
                           vertex_polar: Optional[np.ndarray],
                           out_html: str) -> None:
    traces = []

    (x_pol, y_pol, z_pol,
     x_non, y_non, z_non) = faces_to_wireframe_lines_colored(
        verts_sasa, faces_sasa, vertex_polar
    )

    if x_non is not None and x_non.size > 0:
        traces.append(go.Scatter3d(
            x=x_non, y=y_non, z=z_non,
            mode="lines",
            line=dict(width=WIRE_LINE_WIDTH, color="rgb(0,255,0)"),
            name="SASA wire (non-polar)",
        ))

    if x_pol is not None and x_pol.size > 0:
        traces.append(go.Scatter3d(
            x=x_pol, y=y_pol, z=z_pol,
            mode="lines",
            line=dict(width=WIRE_LINE_WIDTH, color="rgb(255,0,0)"),
            name="SASA wire (polar)",
        ))

    # 原子点
    traces.append(go.Scatter3d(
        x=coords[:, 0],
        y=coords[:, 1],
        z=coords[:, 2],
        mode="markers",
        marker=dict(
            size=5,
            color=[element_color(s) for s in symbols],
            line=dict(width=1, color="black"),
        ),
        name="atoms",
    ))

    fig = go.Figure(
        data=traces,
        layout=go.Layout(
            title=f"{comment} (SASA wireframe with polar segments)",
            scene=dict(aspectmode="data"),
            showlegend=True,
        ),
    )
    fig.write_html(out_html, include_plotlyjs="cdn")
