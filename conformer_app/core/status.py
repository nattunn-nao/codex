from __future__ import annotations
from dataclasses import dataclass
from enum import Enum


class StageState(str, Enum):
    WAIT = "WAIT"
    RUN = "RUN"
    DONE = "DONE"
    ERROR = "ERROR"
    ABORTED = "ABORTED"


@dataclass
class JobProgressRow:
    """
    1 ID 分の進捗テーブル行
    """
    id: str
    smiles: str

    target_confs: int = 0
    generated_confs: int = 0
    mm_optimized_confs: int = 0
    mm_candidates: int = 0
    mm_selected: int = 0

    mm_status: StageState = StageState.WAIT
    uma_status: StageState = StageState.WAIT
    qm_status: StageState = StageState.WAIT

    error_message: str = ""
