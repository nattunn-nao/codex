from __future__ import annotations

from typing import Iterable, List

from conformer_app.core.status import JobProgressRow


def filter_pending_jobs(rows: Iterable[JobProgressRow]) -> List[JobProgressRow]:
    pending: List[JobProgressRow] = []
    for row in rows:
        if (
            row.mm_status.lower() == "done"
            and row.uma_status.lower() == "done"
            and row.qm_status.lower() == "done"
        ):
            continue
        pending.append(row)
    return pending
