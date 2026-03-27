from __future__ import annotations

import argparse
from pathlib import Path

from .config import load_default_config
from .registry import load_surface_registry
from .surface_loader import load_surface
from .workflow import run_batch


def main() -> None:
    parser = argparse.ArgumentParser(description="CMP adsorption workflow CLI")
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--surface", type=str, required=True)
    args = parser.parse_args()

    root = args.root
    cfg = load_default_config(root)
    registry = {e.sid: e for e in load_surface_registry(root)}
    surface_entry = registry[args.surface]
    surface = load_surface(surface_entry.sid, surface_entry.path)

    summaries = run_batch(surface=surface, cfg=cfg, root=root)
    for s in summaries:
        print(f"{s.adsorbate_file}: best={s.best_energy:.4f} eV, report={s.report_csv}")


if __name__ == "__main__":
    main()
