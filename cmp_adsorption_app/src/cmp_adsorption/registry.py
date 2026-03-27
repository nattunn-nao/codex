from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass
class SurfaceEntry:
    sid: str
    label: str
    path: Path


def load_surface_registry(root: Path) -> list[SurfaceEntry]:
    reg_path = root / "configs" / "surface_registry.yaml"
    with reg_path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    out: list[SurfaceEntry] = []
    for item in data.get("surfaces", []):
        out.append(SurfaceEntry(sid=item["id"], label=item.get("label", item["id"]), path=root / item["path"]))
    return out
