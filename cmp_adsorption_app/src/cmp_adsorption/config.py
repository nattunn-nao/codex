from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass
class AppConfig:
    data: dict[str, Any]

    @property
    def input_adsorbates_dir(self) -> Path:
        return Path(self.data.get("input_adsorbates_dir", "adsorbates"))

    @property
    def results_dir(self) -> Path:
        return Path(self.data.get("results_dir", "results"))


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def load_default_config(root: Path) -> AppConfig:
    cfg = load_yaml(root / "configs" / "default_job.yaml")
    return AppConfig(cfg)
