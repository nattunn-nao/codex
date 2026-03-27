from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha1

from ase import Atoms


@dataclass
class RelaxResult:
    structure: Atoms
    energy: float
    tag: str


def _mock_energy(structure: Atoms) -> float:
    payload = structure.get_positions().round(3).tobytes()
    seed = int(sha1(payload).hexdigest()[:8], 16)
    return -50.0 + (seed % 10_000) / 1000.0


def run_gfnff(structures: list[Atoms], use_mock: bool = True) -> list[RelaxResult]:
    results: list[RelaxResult] = []
    for i, st in enumerate(structures):
        if use_mock:
            e = _mock_energy(st)
        else:
            raise NotImplementedError("xTB実行連携は環境依存のため本版ではmockのみ対応")
        results.append(RelaxResult(structure=st, energy=e, tag=f"gfnff_{i:04d}"))
    return results


def prune_by_energy(results: list[RelaxResult], keep_top_n: int) -> list[RelaxResult]:
    return sorted(results, key=lambda x: x.energy)[:keep_top_n]
