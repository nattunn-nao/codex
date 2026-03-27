from __future__ import annotations

from .uma_runner import UmaResult


def rank_structures(results: list[UmaResult]) -> list[UmaResult]:
    return sorted(results, key=lambda r: r.adsorption_energy)
