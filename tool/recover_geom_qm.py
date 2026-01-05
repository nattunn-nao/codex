from __future__ import annotations

import sys
from pathlib import Path

from conformer_app.core.runner import run_geom_qm_extract


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python tool/recover_geom_qm.py <path_to_best_xyz>")
    best_xyz = Path(sys.argv[1])
    if not best_xyz.exists():
        raise SystemExit(f"File not found: {best_xyz}")
    mol_id = best_xyz.stem.split("_")[0]
    qm_dir = best_xyz.parent
    summary_path = run_geom_qm_extract(qm_dir)
    if summary_path:
        target = best_xyz.parent / "Geom_QM_summary.xlsx"
        summary_path.replace(target)
        print(f"Recovered Geom_QM_summary.xlsx for {mol_id}: {target}")


if __name__ == "__main__":
    main()
