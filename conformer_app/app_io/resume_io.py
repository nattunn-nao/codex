# conformer_app/app_io/resume_io.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict
import pandas as pd


@dataclass
class StageFlags:
    mm_status: str = "WAIT"
    uma_status: str = "WAIT"
    qm_status: str = "WAIT"

    @property
    def is_mm_done(self) -> bool:
        return self.mm_status == "DONE"

    @property
    def is_uma_done(self) -> bool:
        return self.uma_status == "DONE"

    @property
    def is_qm_done(self) -> bool:
        return self.qm_status == "DONE"


def load_resume_flags(path_or_file) -> Dict[str, StageFlags]:
    """
    resume_flag.xlsx / progress.xlsx から ID→StageFlags の dict を作る。
    想定カラム:
      - ID
      - MM_status
      - UMA_status
      - QM_status
    """
    df = pd.read_excel(path_or_file)
    cols = {c.lower(): c for c in df.columns}

    id_col = cols.get("id")
    mm_col = cols.get("mm_status")
    uma_col = cols.get("uma_status")
    qm_col = cols.get("qm_status")

    if id_col is None:
        # ID が無いなら再開情報としては使えない
        return {}

    flags: Dict[str, StageFlags] = {}
    for _, row in df.iterrows():
        mol_id = str(row[id_col]).strip()
        if not mol_id:
            continue
        flags[mol_id] = StageFlags(
            mm_status=str(row.get(mm_col, "WAIT")) if mm_col else "WAIT",
            uma_status=str(row.get(uma_col, "WAIT")) if uma_col else "WAIT",
            qm_status=str(row.get(qm_col, "WAIT")) if qm_col else "WAIT",
        )
    return flags
