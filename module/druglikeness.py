'''
Author: Rui Qin
Date: 2025-06-13 20:28:48
LastEditTime: 2025-06-18 20:30:16
Description: 
'''
import importlib.util
from pathlib import Path
import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import QED, Crippen, Mol, rdMolDescriptors, Lipinski
from module.structural import get_rotatable
from utils.preprocess import standard_mol

# Importing the SA_Score module
module_path = Path(RDConfig.RDContribDir) / 'SA_Score' / 'sascorer.py'
spec = importlib.util.spec_from_file_location("sascorer", module_path)
sa = importlib.util.module_from_spec(spec) # type: ignore
spec.loader.exec_module(sa)  # type: ignore


def clogP(mol: Mol) -> float:
    return Crippen.MolLogP(mol) # type: ignore

def hdonors(mol: Mol) -> int:
    return Lipinski.NumHDonors(mol) # type: ignore

def hacceptors(mol: Mol) -> int:
    return Lipinski.NumHAcceptors(mol) # type: ignore

def molwt(mol: Mol) -> float:
    return rdMolDescriptors.CalcExactMolWt(mol)

def qed(mol: Mol) -> float:
    return QED.qed(mol)

def sa_score(mol: Mol) -> float:
    return sa.calculateScore(mol)

def tpsa(mol: Mol) -> float:
    return rdMolDescriptors.CalcTPSA(mol)

def obey_lipinski(mol: Mol) -> int:
    rule_1 = molwt(mol) < 500
    rule_2 = hdonors(mol) <= 5
    rule_3 = hacceptors(mol) <= 10
    rule_4 = (logp := clogP(mol) >= -2) & (logp <= 5)
    rule_5 = get_rotatable(mol) <= 10
    return sum([rule_1, rule_2, rule_3, rule_4, rule_5])

def le_heavyatom(mol: Mol, score: float) -> float:
    """Calculate the ligand efficiency based on the number of heavy atoms and docking score."""
    return mol.GetNumHeavyAtoms() / abs(score) if score <0 else np.nan

def le_mw(mol: Mol, score: float) -> float:
    """Calculate the ligand efficiency based on molecular weight and docking score."""
    return molwt(mol) / abs(score) if score < 0 else np.nan

def calc_le(args:tuple[Mol, float]) -> dict:
    """Calculate the ligand effciency for a given molecule and its docking score.

    Args:
        args (tuple[Chem.Mol, float]): Molecule and its docking score.
    """
    pose, score = args
    mol = standard_mol(pose)
    return {
        'LE_heavyatom': le_heavyatom(mol, score), 
        'LE_mw': le_mw(mol, score)
        }