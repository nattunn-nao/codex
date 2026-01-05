# mm_layer/export.py

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List

from rdkit import Chem

from ase import Atoms

from conformer_app.core.config import OutputConfig
from conformer_app.app_io.paths import mm_dir
from conformer_app.app_io.xyz_io import write_atoms_to_xyz


def _rdkit_conf_to_ase_atoms(
    mol: Chem.Mol,
    conf_id: int,
) -> Atoms:
    """RDKit Mol + conformer ID から ASE Atoms を生成する。"""
    conf = mol.GetConformer(conf_id)

    symbols = []
    positions = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        symbols.append(atom.GetSymbol())
        positions.append((pos.x, pos.y, pos.z))

    atoms = Atoms(symbols=symbols, positions=positions)
    return atoms


def export_selected_conformers(
    mol: Chem.Mol,
    mol_id: str,
    selected_conf_ids: Iterable[int],
    energy_map: Dict[int, float],
    output_config: OutputConfig,
) -> List[str]:
    """
    選抜された conformer を .xyz として出力し、そのファイルパス一覧を返す。

    ファイル名:
        <base_dir>/<ID>/MM/conf001.xyz, conf002.xyz, ...
    コメント行:
        "ID=<ID> conf=<n> confId=<conf_id> E=<energy> kcal/mol"
    """
    base_dir = output_config.base_dir
    out_dir: Path = mm_dir(base_dir, mol_id)

    paths: List[str] = []

    for i, conf_id in enumerate(selected_conf_ids, start=1):
        atoms = _rdkit_conf_to_ase_atoms(mol, conf_id)
        energy = float(energy_map[conf_id])
        comment = (
            f"ID={mol_id} conf={i} confId={conf_id} "
            f"E={energy:.6f} kcal/mol"
        )

        filename = f"conf{i:03d}.xyz"
        path = out_dir / filename

        write_atoms_to_xyz(
            atoms,
            path,
            comment=comment,
            move_to_com=output_config.move_to_com_origin,
        )

        paths.append(str(path))

    return paths
