from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from ase.io import write

from .uma_runner import UmaResult


def ensure_results_dirs(base: Path) -> dict[str, Path]:
    keys = ["inputs", "initial_structures", "gfnff_relaxed", "uma_relaxed", "reports"]
    out = {}
    for k in keys:
        p = base / k
        p.mkdir(parents=True, exist_ok=True)
        out[k] = p
    return out


def save_ranked_results(results: list[UmaResult], report_dir: Path, prefix: str, top_k: int = 5) -> tuple[Path, Path]:
    rows = [{"rank": i + 1, "source_tag": r.source_tag, "adsorption_energy": r.adsorption_energy} for i, r in enumerate(results)]
    df = pd.DataFrame(rows)
    csv_path = report_dir / f"{prefix}_ranking.csv"
    json_path = report_dir / f"{prefix}_summary.json"
    df.to_csv(csv_path, index=False)
    json_path.write_text(json.dumps(rows[:top_k], ensure_ascii=False, indent=2), encoding="utf-8")
    return csv_path, json_path


def dump_top_structures(results: list[UmaResult], out_dir: Path, prefix: str, top_k: int = 5) -> list[Path]:
    paths = []
    for i, item in enumerate(results[:top_k], start=1):
        p = out_dir / f"{prefix}_rank{i:02d}.xyz"
        write(str(p), item.structure)
        paths.append(p)
    return paths
