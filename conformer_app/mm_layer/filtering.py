# mm_layer/filtering.py

from __future__ import annotations

from typing import Dict, Iterable, List

from rdkit import Chem
from rdkit.Chem import AllChem

from conformer_app.core.config import MMConfig


def filter_by_energy_and_window(
    conf_ids: Iterable[int],
    energy_map: Dict[int, float],
    config: MMConfig,
) -> List[int]:
    """
    エネルギー昇順にソートし、
    - ΔE ≤ equal_de_threshold_kcal を同値として 1 つに代表化
    - Emin + energy_window_kcal を越えたら打ち切り
    という条件で候補 conformer のリストを返す。
    """
    conf_ids = list(conf_ids)
    if not conf_ids:
        return []

    # エネルギーでソート
    sorted_conf_ids = sorted(conf_ids, key=lambda cid: energy_map[cid])
    Emin = energy_map[sorted_conf_ids[0]]

    equal_de = float(config.equal_de_threshold_kcal)
    window = float(config.energy_window_kcal)

    accepted: List[int] = []

    for cid in sorted_conf_ids:
        e = energy_map[cid]

        # エネルギー窓 (E <= Emin + window)
        if e > Emin + window:
            break

        # 同値判定: すでに受理した構造との ΔE が equal_de 以下ならスキップ
        is_equivalent = any(
            abs(e - energy_map[prev_cid]) <= equal_de for prev_cid in accepted
        )
        if is_equivalent:
            continue

        accepted.append(cid)

    return accepted


def select_by_rmsd(
    mol: Chem.Mol,
    candidate_conf_ids: List[int],
    config: MMConfig,
) -> List[int]:
    """
    剛体 Align & RMSD による多様性選抜。

    - 最もエネルギーが低い構造（candidate_conf_ids[0] 前提）を参照構造とする
    - 参照との差 RMSD >= rmsd_min を満たすものを、
      エネルギー順（candidate_conf_ids の順と同じ）に追加していく
    - 最大 max_selected_confs 個まで
    """
    if not candidate_conf_ids:
        return []

    rmsd_min = float(config.rmsd_min)
    max_selected = int(config.max_selected_confs)

    # candidate_conf_ids は既にエネルギー昇順と想定
    ref_id = candidate_conf_ids[0]
    selected: List[int] = [ref_id]

    for cid in candidate_conf_ids[1:]:
        if len(selected) >= max_selected:
            break

        rms = AllChem.GetBestRMS(mol, mol, ref_id, cid)
        if rms >= rmsd_min:
            selected.append(cid)

    return selected
