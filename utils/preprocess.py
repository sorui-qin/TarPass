'''
Author: Rui Qin
Date: 2025-03-16 15:03:08
LastEditTime: 2025-04-14 20:47:56
Description: 
'''
from typing import Tuple, List
from collections import defaultdict
from pathlib import Path
from utils.logger import project_logger
from utils.io import read_sdf, read_smi
from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') # type: ignore

def to_mols(smiles:list) -> List[Chem.Mol]:
    return [Chem.MolFromSmiles(smile) for smile in smiles]

def to_smiles(mols:list) -> List[str]:
    return [Chem.MolToSmiles(mol) for mol in mols]

def standard_mol(mol:Chem.Mol) -> Chem.Mol:
    """Reset the molecule to *standard* form without conformation.  
    Not same as `Chem.MolStandardize.rdMolStandardize.Cleanup((Mol)`
    """
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))

def _sanitize_valid(mol:Chem.Mol, idx:int) -> bool:
    mol.SetProp('_Name', f'Mol ID {idx}') # set the name of the molecule
    try:
        Chem.SanitizeMol(mol)
        return True
    except:
        return False

def _deduplicate_3D(mols, indices):
    kept = []
    for elem in indices:
        conflict = False
        for kept_elem in kept:
            if CalcRMS(mols[elem], mols[kept_elem]) == 0:
                conflict = True
                break
        if not conflict:
            kept.append(elem)
    return kept


class Preprocess():
    def __init__(self, mols):
        self.mols = mols
        self.mols_num = len(mols)
        
    def valid(self) -> Tuple[List[str], List[Chem.Mol]]:
        valids = [mol for idx, mol in enumerate(self.mols) if _sanitize_valid(mol, idx)]
        project_logger.info(f'Valid molecules: {len(valids)} out of {self.mols_num}')
        return to_smiles(valids), valids
    
    def unique(self, if_valid=True) -> Tuple[List[str], List[Chem.Mol]]:
        """Check the uniqueness of the molecule list.   
        It will return a deduplicated list of canonical SMILES `unique_smis` and their corresponding original Mol object lists `unique_mols`.  
        It is important to note that if there are different 3D conformations with the same canonical SMILES,
        the Mol objects will be assessed based on RMSD to determine if they represent the same conformation.   
        If no duplicates are found, all conformations will be retained.

        Args:
            if_valid (bool, optional): Using molecules list passed the valid check. Defaults to True.
        """
        mols = self.mols if not if_valid else self.valid()[1]
        unique_di = defaultdict(list)
        for i, mol in enumerate(mols):
            smi = Chem.MolToSmiles(mol)
            unique_di[smi].append(i)

        consider_3D = False
        unique_smis, unique_mols = [], []
        for smi, indices in unique_di.items():
            unique_smis.append(smi)
            if mols[indices[0]].GetNumConformers(): # Consider conformer uniqueness # type: ignore
                consider_3D = True
                if len(indices) > 1:
                    indices = _deduplicate_3D(mols, indices)
            else: # If no conformer, only keep the first one
                indices = [indices[0]]
            unique_mols.extend([mols[idx] for idx in indices])
        
        project_logger.info(f'Unique SMILES: {len(unique_smis)} out of {self.mols_num}')
        if consider_3D:
            project_logger.info(f'3D Conformation checked...')
            project_logger.info(f'Unique Conformation: {len(unique_mols)} out of {self.mols_num}')
        return unique_smis, unique_mols


def read_in(target_dir, unique=True) -> Tuple[List[str], List[Chem.Mol]]:
    """Read in molecules from the target directory.
    Args:
        target_dir (Path): Path to the target directory.
        unique (bool, optional): Running uniqueness process.

    Returns:
        Tuple[List[str], List[Chem.Mol]]: Processed SMILES list `smis` and Mol list `mols`.
    """
    read_dir = Path(target_dir) #/'generated' # Modify here when the read-in path is different

    sdf_files = list(read_dir.glob('*.sdf'))
    if sdf_files: #SDF files will be prioritized for reading
        project_logger.info(f"Found {len(sdf_files)} SDF file(s). Reading SDF files...")
        read_files, reader = sdf_files, read_sdf
    else: # If no SDF files found, read SMILES
        project_logger.info("No SDF files found. Reading SMILES instead...")
        read_files, reader = list(Path(read_dir).iterdir()), read_smi

    mols = []
    for f in read_files:
        if f.is_dir():
            continue
        try:
            readin = reader(f)
            mols += readin if isinstance(readin, list) else [readin]
        except Exception as e:
            project_logger.warning(f"Error reading {f}: {e}")
    process = Preprocess(mols)
    return process.unique() if unique else process.valid()