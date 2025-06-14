'''
Author: Rui Qin
Date: 2025-06-13 20:28:48
LastEditTime: 2025-06-14 14:52:25
Description: 
'''
import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from utils.preprocess import standard_mol


def molwt(mol: Chem.Mol) -> float:
    """Calculate the molecular weight of a given molecule."""
    return CalcExactMolWt(mol)

def calc_le(args:tuple[Chem.Mol, float]) -> dict:
    pose, score = args
    mol = standard_mol(pose)
    if score <= 0:
        le_heavyatom = mol.GetNumHeavyAtoms() / abs(score) if abs(score) > 0 else np.nan
        le_mw = molwt(mol) / abs(score) if abs(score) > 0 else np.nan
    return {'LE_heavyatom': le_heavyatom, 'LE_mw': le_mw}