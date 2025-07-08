'''
Author: Rui Qin
Date: 2025-03-16 15:03:08
LastEditTime: 2025-07-08 16:21:10
Description: 
'''
from copy import deepcopy
from pathlib import Path
from typing import DefaultDict, Iterable, Literal
from collections import defaultdict
from itertools import islice
from utils.logger import project_logger
from utils.io import read_sdf, read_strings
from utils.constant import DASHLINE
from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') # type: ignore

def to_mols(smiles:Iterable) -> list[Chem.Mol]:
    return [mol for smile in smiles if (mol:=Chem.MolFromSmiles(smile))]

def to_smiles(mols:Iterable) -> list[str]:
    return [smi for mol in mols if (smi:=Chem.MolToSmiles(mol, kekuleSmiles=True))]

def to_inchi(mols:Iterable) -> list:
    return [Chem.MolToInchi(mol) for mol in mols if mol]

def standard_mol(mol:Chem.Mol) -> Chem.Mol:
    """Reset the molecule to *standard* form without conformation.  
    Not same as `Chem.MolStandardize.rdMolStandardize.Cleanup((Mol)`
    """
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))

def canonical_smiles(smiles:str) -> str:
    """Reset the SMILES to canonical form.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol else smiles
    except:
        return smiles

def sanitize_valid(mol:Chem.Mol, idx:int) -> bool:
    """Check if a mol object is valid.
    """
    if mol:
        mol.SetProp('_Name', f'{idx}') # set the name of the molecule
        if not Chem.MolToSmiles(mol): # Check the empty Mol object.
            return False
        try:
            molcopy = deepcopy(mol)
            Chem.SanitizeMol(molcopy)
            return True
        except:
            return False
    else:
        return False
    
def smiles_valid(smi:str, idx:int) -> Chem.Mol|None:
    """Check if a SMILES string is valid.
    """
    mol = Chem.MolFromSmiles(smi)
    mol.SetProp('_Name', f'{idx}') # set the name of the molecule
    if mol is None:
        return None
    return mol

def conformation_check(mols:list[Chem.Mol]) -> bool:
    """Check if all molecules have 3D conformations.
    Returns:
        bool: True if all molecules have 3D conformations, False otherwise.
    """
    confs = [1 if mol.GetNumConformers() else 0 for mol in mols]
    if any(confs) + any(confs) == 0:
        return False
    elif any(confs) + any(confs) == 1:
        raise RuntimeError('Some molecules have not 3D conformation. Please check the input molecules.')
    return True

def check_duplicate3D(mols:list[Chem.Mol]) -> list[Chem.Mol]:
    """Check if the new molecule is a duplicate of any existing 3D conformations.
    """
    if not conformation_check(mols):
        return mols  # If no 3D conformations, return the original list
    deplicate = []
    for mol in mols:
        is_duplicate = False
        for unique_mol in deplicate:
            if CalcRMS(mol, unique_mol) == 0:
                is_duplicate = True
                break
        if not is_duplicate:
            deplicate.append(mol)
    return deplicate


class Preprocess():
    def __init__(self, readins:list, format:Literal['sdf', 'smi']):
        self.readins = readins
        self.num = len(readins)
        self.format = format
    
    def valid(self) -> list[Chem.Mol]:
        if self.format == 'sdf':
            valids = [mol for idx, mol in enumerate(self.readins) if sanitize_valid(mol, idx)]
        elif self.format == 'smi':
            valids = [mol for idx, smi in enumerate(self.readins) if (mol:=smiles_valid(smi, idx))]
        project_logger.info(f'Valid molecules: {len(valids)} out of {self.num}')
        return valids

    def unique(self) -> DefaultDict[str, list[Chem.Mol]]:
        """Check the uniqueness of the molecule list.
        Returns:
            DefaultDict[str, List[Chem.Mol]]: A dictionary with unique SMILES as keys and a list of corresponding Mol objects as values.
        """
        mols = self.valid()
        unique_di = defaultdict(list)
        for mol in mols:
            smi = Chem.MolToSmiles(mol, kekuleSmiles=True)
            unique_di[smi].append(mol)
        project_logger.info(f'Unique SMILES: {len(unique_di.items())} out of {len(mols)}')
        return unique_di

def read_in(target_dir, num_thres=1000, isomers=False) -> tuple[list[str], list[Chem.Mol]]:
    """Read in molecules from the target directory. Return the duplicate SMILES list and Mol list.
    Args:
        target_dir (Path): Path to the target directory.
        num_thres (int, optional): Threshold for the number of molecules. Defaults to 1000. If none, all molecules will be read in.
        isomers (bool, optional): Whether to consider isomers. If True, the Mol list will contain isomers, which will make the list slightly longer than num_thres. Defaults to False.  
    Returns:
        Tuple[List[str], List[Chem.Mol]]: Processed SMILES list `smis` and Mol list `mols`.
    """
    read_dir = Path(target_dir) #/'generated' # Modify here when the read-in path is different

    # Check all readable files
    sdf_files = sorted(read_dir.glob('*.sdf'))
    if sdf_files: #SDF files will be prioritized for reading
        project_logger.info(f"Found {len(sdf_files)} SDF file(s). Reading SDF files...")
        read_files, reader, format = sdf_files, read_sdf, 'sdf'
    else: # If no SDF files found, read SMILES
        project_logger.info("No SDF files found. Reading SMILES instead...")
        read_files, reader, format = sorted(read_dir.iterdir()), read_strings, 'smi'

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
    if num_thres is None or num_thres <= 0:
        num_thres = len(process.readins)
    # Limit the number of molecules to num_thres
    processed_li = list(islice(process.unique().items(), num_thres))
    proc_smis, proc_mols = [], []

    # Reading SMILES or not considering isomers
    if format == 'smi' or not isomers:
        for smi, mols in processed_li:
            proc_smis.append(smi)
            proc_mols.append(mols[0])
    # Considering isomers:
    if isomers and format == 'sdf':
        project_logger.info(f'Checking 3D Conformation...')
        for smi, mols in processed_li:
            proc_smis.append(smi)
            proc_mols.extend(check_duplicate3D(mols))
        project_logger.info(f'3D Conformation checked, {len(proc_mols)} unique 3D conformations found.')

    if len(proc_smis) < num_thres:
        project_logger.warning(f"Not enough unique molecules found. Read with {len(proc_smis)} unique molecules, but expected {num_thres}.")
    else:
        project_logger.info(f"Read in procession is done, top {num_thres} unique molecules are selected.")
    project_logger.info(DASHLINE)
    return proc_smis, proc_mols