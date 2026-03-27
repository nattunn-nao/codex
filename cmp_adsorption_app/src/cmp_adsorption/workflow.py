from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ase.io import write

from .adsorbate_loader import list_adsorbate_xyz, load_adsorbate
from .config import AppConfig
from .geometry_filter import geometry_pass
from .gfnff_runner import prune_by_energy, run_gfnff
from .pose_generator import PoseConfig, generate_initial_poses
from .ranker import rank_structures
from .reporter import dump_top_structures, ensure_results_dirs, save_ranked_results
from .surface_loader import SurfaceModel
from .uma_runner import run_uma_oc25


@dataclass
class WorkflowSummary:
    adsorbate_file: str
    n_initial: int
    n_after_geom: int
    n_after_gfnff: int
    best_energy: float
    report_csv: Path


def run_single_adsorbate(surface: SurfaceModel, ads_path: Path, cfg: AppConfig, root: Path) -> WorkflowSummary:
    params = cfg.data
    dirs = ensure_results_dirs(root / cfg.results_dir)

    ads = load_adsorbate(ads_path)
    write(str(dirs["inputs"] / ads_path.name), ads)

    pose_cfg = PoseConfig(**params["pose"])
    poses = generate_initial_poses(surface.slab, ads, pose_cfg)

    geom_ok = [p for p in poses if geometry_pass(p, len(surface.slab), params["filters"]["min_distance"])]

    gfnff = run_gfnff(geom_ok, use_mock=params.get("mock_gfnff", True))
    kept = prune_by_energy(gfnff, keep_top_n=params["pruning"]["keep_top_n_after_gfnff"])

    uma_input = [(r.tag, r.structure) for r in kept]
    uma = run_uma_oc25(uma_input, use_mock=params.get("mock_uma", True))
    ranked = rank_structures(uma)

    prefix = f"{surface.sid}_{ads_path.stem}"
    csv_path, _ = save_ranked_results(ranked, dirs["reports"], prefix=prefix, top_k=params["ranking"]["top_k"])
    dump_top_structures(ranked, dirs["uma_relaxed"], prefix=prefix, top_k=params["ranking"]["top_k"])

    for i, st in enumerate(poses[: min(20, len(poses))]):
        write(str(dirs["initial_structures"] / f"{prefix}_init_{i:03d}.xyz"), st)
    for i, st in enumerate([k.structure for k in kept]):
        write(str(dirs["gfnff_relaxed"] / f"{prefix}_gfnff_{i:03d}.xyz"), st)

    return WorkflowSummary(
        adsorbate_file=ads_path.name,
        n_initial=len(poses),
        n_after_geom=len(geom_ok),
        n_after_gfnff=len(kept),
        best_energy=ranked[0].adsorption_energy if ranked else float("nan"),
        report_csv=csv_path,
    )


def run_batch(surface: SurfaceModel, cfg: AppConfig, root: Path) -> list[WorkflowSummary]:
    ads_files = list_adsorbate_xyz(root / cfg.input_adsorbates_dir)
    summaries = []
    for ads_path in ads_files:
        summaries.append(run_single_adsorbate(surface=surface, ads_path=ads_path, cfg=cfg, root=root))
    return summaries
