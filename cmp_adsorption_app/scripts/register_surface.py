#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path.cwd())
    p.add_argument("--id", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--path", required=True, help="Path relative to root")
    args = p.parse_args()

    registry_path = args.root / "configs" / "surface_registry.yaml"
    data = yaml.safe_load(registry_path.read_text(encoding="utf-8")) or {"surfaces": []}
    data.setdefault("surfaces", []).append({"id": args.id, "label": args.label, "path": args.path})
    registry_path.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(f"registered: {args.id}")


if __name__ == "__main__":
    main()
