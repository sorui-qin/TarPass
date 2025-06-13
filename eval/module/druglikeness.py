from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt

def molwt(mol: Chem.Mol) -> float:
    """Calculate the molecular weight of a given molecule."""
    return CalcExactMolWt(mol)