from __future__ import annotations

import os

# OpenMP ランタイムの二重ロード対策（fairchem/torch より前に設定）
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from fairchem.core import pretrained_mlip, FAIRChemCalculator  # type: ignore

from conformer_app.core.config import UMAConfig


def create_uma_calculator(cfg: UMAConfig) -> FAIRChemCalculator:
    """
    UMAConfig から fairchem の Calculator を構築する。
    """
    model_name = cfg.model_name
    device = cfg.device
    task_name = cfg.task_name or "omol"

    predictor = pretrained_mlip.get_predict_unit(model_name, device=device)
    calc = FAIRChemCalculator(predictor, task_name=task_name)
    return calc
