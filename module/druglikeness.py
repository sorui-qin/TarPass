'''
Author: Rui Qin
Date: 2025-06-13 20:28:48
LastEditTime: 2025-06-27 16:40:38
Description: 
'''
import importlib.util
from pathlib import Path
from typing import Union
import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import QED, Crippen, Mol, rdMolDescriptors
from utils.preprocess import standard_mol
from module.filter import match_pains, match_SureChEMBL, match_glaxo

# Importing the SA_Score module
module_path = Path(RDConfig.RDContribDir) / 'SA_Score' / 'sascorer.py'
spec = importlib.util.spec_from_file_location("sascorer", module_path)
sa = importlib.util.module_from_spec(spec) # type: ignore
spec.loader.exec_module(sa)  # type: ignore


class DruglikenessCalculator:
    def __init__(self):
        self.properties = {}
    
    def calculate_all(self, mol:Chem.Mol) -> dict[str, Union[float, int]]:
        if mol is None:
            return {}
        
        self.properties = {
            # Physicochemical Properties
            'clogP': Crippen.MolLogP(mol), # type: ignore
            'hdonors': rdMolDescriptors.CalcNumHBD(mol),
            'hacceptors': rdMolDescriptors.CalcNumHBA(mol),
            'molwt': rdMolDescriptors.CalcExactMolWt(mol), 
            'qed': QED.qed(mol),
            'tpsa': rdMolDescriptors.CalcTPSA(mol),
            'sa_score': sa.calculateScore(mol),
            
            # Global Structural Descriptors
            'heavy_atoms': mol.GetNumHeavyAtoms(),
            'heteroatoms': rdMolDescriptors.CalcNumHeteroatoms(mol),
            'all_common': common_atoms(mol),
            'csp3': rdMolDescriptors.CalcFractionCSP3(mol),
            'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
            'spiro_atoms': rdMolDescriptors.CalcNumSpiroAtoms(mol),
            'chiral_atoms': len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
            'completeness': len(Chem.GetMolFrags(mol)) == 1,
        }
        
        # Derived Properties
        self.properties.update({
            'heteroatom_ratio': self.heteroatom_ratio(),
            'flexibility': self.flexibility(),
            'lipinski': self.lipinski(),
        })

        # Structural Alerts
        self.properties.update({
            'PAINS_alert': match_pains(mol),
            'SureChEMBL_alert': match_SureChEMBL(mol),
            'Glaxo_alert': match_glaxo(mol),
        })
        return self.properties
    
    def heteroatom_ratio(self) -> float:
        heavy_atoms = self.properties.get('heavy_atoms', 0)
        return (self.properties.get('heteroatoms', 0) / heavy_atoms 
                if heavy_atoms > 0 else np.nan)
    
    def flexibility(self) -> float:
        total_bonds = self.properties.get('heavy_atoms', 1) - 1  # 近似总键数
        rot_bonds = self.properties.get('rotatable_bonds', 0)
        rigid_bonds = max(1, total_bonds - rot_bonds)
        return rot_bonds / rigid_bonds
    
    def lipinski(self) -> int:
        rule_1 = self.properties.get('molwt', 0) >= 500
        rule_2 = self.properties.get('hdonors', 0) <=5
        rule_3 = self.properties.get('hacceptors', 0) <= 10
        logp = self.properties.get('clogP', 0)
        rule_4 = (logp >= -2) & (logp <= 5)
        rule_5 = self.properties.get('rotatable_bonds', 0) <= 101
        return sum([rule_1, rule_2, rule_3, rule_4, rule_5])
    
#### Common Atoms #####

def common_atoms(mol:Mol) -> bool:
    """Check if the molecule doesn't contain uncommon atoms.  
    Namely, the atom does not belong to the following categories:  
    core elements `C, H, O, N`, halogens `F, Cl, Br, I`, or common functional group elements `P, S, B`.
    """
    common = {'C', 'H', 'O', 'N', 'F', 'Cl', 'Br', 'I', 'P', 'S', 'B'}
    return all(atom.GetSymbol() in common for atom in mol.GetAtoms())

#### Ligand Effciency #####

def le_heavyatom(mol:Mol, score: float) -> float:
    """Calculate the ligand efficiency based on the number of heavy atoms and docking score."""
    return abs(score) / mol.GetNumHeavyAtoms() if score < 0 else np.nan

def le_mw(mol:Mol, score: float) -> float:
    """Calculate the ligand efficiency based on molecular weight and docking score."""
    return abs(score) / rdMolDescriptors.CalcExactMolWt(mol) if score < 0 else np.nan

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