# Adapted from PoseBusters: https://github.com/maabuu/posebusters
"""Module to check intermolecular distances between ligand and protein."""

from __future__ import annotations

from typing import Any, Iterable

import numpy as np
import pandas as pd
from rdkit.Chem.rdchem import Atom, GetPeriodicTable, Mol

_periodic_table = GetPeriodicTable()

def check_intermolecular_distance(  # noqa: PLR0913
    mol_pred: Mol,
    mol_cond: Mol,
    radius_type: str = "vdw",
    radius_scale: float = 1.0,
    clash_cutoff: float = 0.75,
    ignore_types: set[str] = {"hydrogens"},
    max_distance: float = 5.0,
    search_distance: float = 6.0,
) -> dict[str, Any]:
    """Check that predicted molecule is not too close and not too far away from conditioning molecule.

    Args:
        mol_pred: Predicted molecule (docked ligand) with one conformer.
        mol_cond: Conditioning molecule (protein) with one conformer.
        radius_type: Type of atomic radius to use. Possible values are "vdw" (van der Waals) and "covalent".
            Defaults to "vdw".
        radius_scale: Scaling factor for the atomic radii. Defaults to 0.8.
        clash_cutoff: Threshold for how much the atoms may overlap before a clash is reported. Defaults
            to 0.05.
        ignore_types: Which types of atoms to ignore in mol_cond. Possible values to include are "hydrogens", "protein",
            "organic_cofactors", "inorganic_cofactors", "waters". Defaults to {"hydrogens"}.
        max_distance: Maximum distance (in Angstrom) predicted and conditioning molecule may be apart to be considered
            as valid. Defaults to 5.0.

    Returns:
        PoseBusters results dictionary.
    """
    coords_ligand = mol_pred.GetConformer().GetPositions()
    coords_protein = mol_cond.GetConformer().GetPositions()

    atoms_ligand = np.array([a.GetSymbol() for a in mol_pred.GetAtoms()])
    atoms_protein_all = np.array([a.GetSymbol() for a in mol_cond.GetAtoms()])

    #idxs_ligand = np.array([a.GetIdx() for a in mol_pred.GetAtoms()])
    idxs_protein = np.array([a.GetIdx() for a in mol_cond.GetAtoms()])

    mask = [a.GetSymbol() != "H" for a in mol_pred.GetAtoms()]
    coords_ligand = coords_ligand[mask, :]
    atoms_ligand = atoms_ligand[mask]
    #mask_ligand_idxs = idxs_ligand[mask]
    if ignore_types:
        mask = get_atom_type_mask(mol_cond, ignore_types)
        coords_protein = coords_protein[mask, :]
        atoms_protein_all = atoms_protein_all[mask]
        mask_protein_idxs = idxs_protein[mask]

    # get radii
    radius_ligand = _get_radii(atoms_ligand, radius_type)
    radius_protein_all = _get_radii(atoms_protein_all, radius_type)

    # select atoms that are close to ligand to check for clash
    distances_all = _pairwise_distance(coords_ligand, coords_protein)
    mask_protein = distances_all.min(axis=0) <= search_distance
    distances = distances_all[:, mask_protein]
    radius_protein = radius_protein_all[mask_protein]
    #atoms_protein = atoms_protein_all[mask_protein]
    mask_protein_idxs = mask_protein_idxs[mask_protein]

    radius_sum = radius_ligand[:, None] + radius_protein[None, :]
    relative_distance = distances / radius_sum
    violations = relative_distance < 1 / radius_scale

    if distances.size > 0:
        violations[np.unravel_index(distances.argmin(), distances.shape)] = True  # add smallest distances as info
        violations[np.unravel_index(relative_distance.argmin(), relative_distance.shape)] = True
    violation_ligand, violation_protein = np.where(violations)
    #reverse_ligand_idxs = mask_ligand_idxs[violation_ligand]
    #reverse_protein_idxs = mask_protein_idxs[violation_protein]

    # collect details around those violations in a dataframe
    details = pd.DataFrame()
    #details["ligand_atom_id"] = reverse_ligand_idxs
    #details["protein_atom_id"] = reverse_protein_idxs
    #details["ligand_element"] = [atoms_ligand[i] for i in violation_ligand]
    #details["protein_element"] = [atoms_protein[i] for i in violation_protein]
    details["ligand_vdw"] = [radius_ligand[i] for i in violation_ligand]
    details["protein_vdw"] = [radius_protein[i] for i in violation_protein]
    details["sum_radii"] = details["ligand_vdw"] + details["protein_vdw"]
    details["distance"] = distances[violation_ligand, violation_protein]
    details["sum_radii_scaled"] = details["sum_radii"] * radius_scale
    details["relative_distance"] = details["distance"] / details["sum_radii_scaled"]
    details["clash"] = details["relative_distance"] < clash_cutoff

    results = {
        "no_clashes": not details["clash"].any(),
    }

    return results


def _pairwise_distance(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.linalg.norm(x[:, None, :] - y[None, :, :], axis=-1)


def _get_radii(atoms: np.ndarray, radius_type: str) -> np.ndarray:
    if radius_type == "vdw":
        return np.array([_periodic_table.GetRvdw(a) for a in atoms])
    elif radius_type == "covalent":
        return np.array([_periodic_table.GetRcovalent(a) for a in atoms])
    else:
        raise ValueError(f"Unknown radius type {radius_type}. Valid values are 'vdw' and 'covalent'.")
    

def get_atom_type_mask(mol: Mol, ignore_types: Iterable[str]) -> list[bool]:
    """Get mask for atoms to keep."""
    ignore_types = set(ignore_types)
    unsupported = ignore_types - {"hydrogens", "protein", "organic_cofactors", "inorganic_cofactors", "waters"}
    if unsupported:
        raise ValueError(f"Ignore types {unsupported} not supported.")

    if mol.GetAtomWithIdx(0).GetPDBResidueInfo() is None:
        #logger.warning("No PDB information found. Assuming organic molecule.")
        pass

    ignore_h = "hydrogens" in ignore_types
    ignore_protein = "protein" in ignore_types
    ignore_org_cof = "organic_cofactors" in ignore_types
    ignore_inorg_cof = "inorganic_cofactors" in ignore_types
    ignore_water = "waters" in ignore_types

    return [
        _keep_atom(a, ignore_h, ignore_protein, ignore_org_cof, ignore_inorg_cof, ignore_water) for a in mol.GetAtoms()
    ]


_inorganic_cofactor_elements = {
    "Li",
    "Be",
    "Na",
    "Mg",
    "Cl",
    "K",
    "Ca",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Br",
    "Rb",
    "Mo",
    "Cd",
}
_inorganic_cofactor_ccd_codes = {
    "FES",
    "MOS",
    "PO3",
    "PO4",
    "PPK",
    "SO3",
    "SO4",
    "VO4",
}
_water_ccd_codes = {
    "HOH",
}

def _keep_atom(  # noqa: PLR0913, PLR0911
    atom: Atom, ignore_h: bool, ignore_protein: bool, ignore_org_cof: bool, ignore_inorg_cof: bool, ignore_water: bool
) -> bool:
    """Whether to keep atom for given ignore flags."""
    symbol = atom.GetSymbol()
    if ignore_h and symbol == "H":
        return False

    if ignore_inorg_cof and symbol in _inorganic_cofactor_elements:
        return False

    # if loaded from PDB file, we can use the residue names and the hetero flag
    info = atom.GetPDBResidueInfo()
    if info is None:
        if ignore_org_cof:
            return False
        return True

    is_hetero = info.GetIsHeteroAtom()
    if ignore_protein and not is_hetero:
        return False

    residue_name = info.GetResidueName()
    if ignore_water and residue_name in _water_ccd_codes:
        return False

    if ignore_inorg_cof and residue_name in _inorganic_cofactor_ccd_codes:
        return False

    return True

from utils.docking import centriod

def centriod_distance(mol_pred: Mol, mol_ref: Mol):
    """Calculate the distance between the centroids of two molecules."""
    centroid_pred = centriod(mol_pred)
    centroid_cond = centriod(mol_ref)
    return np.linalg.norm(
        np.array(centroid_pred) - np.array(centroid_cond)
        )