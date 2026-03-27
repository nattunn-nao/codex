from __future__ import annotations

from dataclasses import dataclass
from hashlib import md5

from ase import Atoms


@dataclass
class UmaResult:
    structure: Atoms
    adsorption_energy: float
    source_tag: str


def _mock_uma_energy(structure: Atoms) -> float:
    payload = structure.get_positions().round(4).tobytes()
    seed = int(md5(payload).hexdigest()[:8], 16)
    return -5.0 + (seed % 5000) / 1000.0


def run_uma_oc25(candidates: list[tuple[str, Atoms]], use_mock: bool = True) -> list[UmaResult]:
    out: list[UmaResult] = []
    for tag, structure in candidates:
        if use_mock:
            e = _mock_uma_energy(structure)
        else:
            raise NotImplementedError("UMA/OC25推論呼び出しは環境依存のため本版ではmockのみ対応")
        out.append(UmaResult(structure=structure, adsorption_energy=e, source_tag=tag))
    return out
