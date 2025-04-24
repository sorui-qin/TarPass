'''
Author: Rui Qin
Date: 2025-03-16 15:03:08
LastEditTime: 2025-04-24 11:55:26
Description: 
'''
from pathlib import Path
from typing import Tuple, List, DefaultDict
from collections import defaultdict
from itertools import islice
from utils.logger import project_logger
from utils.io import read_sdf, read_smi, read_strings
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

def sanitize_valid(mol:Chem.Mol, idx:int) -> bool:
    """Check if a mol object is valid.
    """
    mol.SetProp('_Name', f'MolID {idx}') # set the name of the molecule
    try:
        Chem.SanitizeMol(mol)
        return True
    except:
        return False
    
def smiles_valid(smi:str, idx:int) -> Chem.Mol|None:
    """Check if a SMILES string is valid.
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol.SetProp('_Name', f'MolID {idx}') # set the name of the molecule
    return mol

def _deduplicate_3D(mols):
    """Deprecated.
    """
    unique = []
    for mol in mols:
        is_duplicate = False
        for unique_mol in unique:
            if CalcRMS(mol, unique_mol) == 0:
                is_duplicate = True
                break
        if not is_duplicate:
            unique.append(mol)
    return unique

def check_duplicate3D(mols:List[Chem.Mol], new_mol:Chem.Mol) -> bool:
    """Check if the new molecule is a duplicate of any existing 3D conformations.
    """
    for mol in mols:
        if CalcRMS(mol, new_mol) == 0:
            return False
    return True


class Preprocess():
    def __init__(self, readins, format):
        self.readins = readins
        self.num = len(readins)
        self.format = format
        
    def valid(self) -> Tuple[List[str], List[Chem.Mol]]:
        if self.format == 'sdf':
            valids = [mol for idx, mol in enumerate(self.readins) if sanitize_valid(mol, idx)]
        elif self.format == 'smi':
            valids = [mol for idx, smi in enumerate(self.readins) if (mol:=smiles_valid(smi, idx))]
        project_logger.info(f'Valid molecules: {len(valids)} out of {self.num}')
        return to_smiles(valids), valids

    def unique(self) -> DefaultDict[str, List[Chem.Mol]]:
        """Check the uniqueness of the molecule list.
        Args:
            if_valid (bool, optional): Using molecules list passed the valid check. Defaults to True.
        """
        check_3D = []
        mols = self.valid()[1]
        unique_di = defaultdict(list)
        for mol in mols:
            conf = 1 if mol.GetNumConformers() else 0
            smi = Chem.MolToSmiles(mol)
            if unique_di[smi]:
                if conf == 0 or (conf == 1 and not check_duplicate3D(unique_di[smi], mol)):
                    continue
            unique_di[smi].append(mol)
            check_3D.append(conf)

        # Check if all molecules are 3D or 2D conformations
        consider_3D = all(check_3D)
        if not consider_3D and sum(check_3D):
            raise RuntimeError('Some molecules are 3D conformations, but not all of them. Please check the input molecules.')
        
        project_logger.info(f'Unique SMILES: {len(unique_di.items())} out of {len(mols)}')
        if consider_3D:
            project_logger.info(f'3D Conformation checked...')
            project_logger.info(f'Unique Conformation: {sum(len(v) for v in unique_di.values())} out of {len(mols)}')
        return unique_di


def read_in(target_dir, num_thres=1000) -> Tuple[List[str], List[Chem.Mol]]:
    """Read in molecules from the target directory. Return the duplicate SMILES list and Mol list.
    Args:
        target_dir (Path): Path to the target directory.
        num_thres (int, optional): Threshold for the number of molecules. Defaults to 1000.

    Returns:
        Tuple[List[str], List[Chem.Mol]]: Processed SMILES list `smis` and Mol list `mols`.
    """
    read_dir = Path(target_dir) #/'generated' # Modify here when the read-in path is different

    # Check all readable files
    sdf_files = list(read_dir.glob('*.sdf'))
    if sdf_files: #SDF files will be prioritized for reading
        project_logger.info(f"Found {len(sdf_files)} SDF file(s). Reading SDF files...")
        read_files, reader, format = sdf_files, read_sdf, 'sdf'
    else: # If no SDF files found, read SMILES
        project_logger.info("No SDF files found. Reading SMILES instead...")
        read_files, reader, format = list(Path(read_dir).iterdir()), read_strings, 'smi'

    # Read in the files
    readins = []
    for f in read_files:
        if f.is_dir():
            continue
        try:
            readin = reader(f)
            readins += readin if isinstance(readin, list) else [readin]
        except Exception as e:
            project_logger.warning(f"Error reading {f}: {e}")

    # Unique and valid process
    process = Preprocess(readins, format)
    processed_di = dict(islice(process.unique().items(), num_thres)) # Limit the number of molecules to num_thres
    proc_smis, proc_mols = [], []
    for smi, mols in processed_di.items():
        proc_smis.append(smi)
        proc_mols.extend(mols)
    project_logger.info(f"Read-in procession is done, top {num_thres} unique molecules are selected.")
    return proc_smis, proc_mols