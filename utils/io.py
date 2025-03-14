'''
Author: Rui Qin
Date: 2025-03-08 15:38:31
LastEditTime: 2025-03-14 21:16:19
Description: 
'''
import os
import pickle
import tempfile
from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS
from typing import Optional, Tuple, List
from collections import defaultdict
from collections.abc import Iterable
from utils.logger import project_logger
from contextlib import contextmanager
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def write_pkl(pkl_file, data):
    with open(pkl_file, 'wb') as fi:
        pickle.dump(data, fi)

def read_pkl(pkl_file):
    data = []
    with open(pkl_file, 'rb') as f:
        while True:
            try:
                aa = pickle.load(f)
                data.extend(aa)
            except EOFError:
                break
    return data

def read_sdf(sdf_file, sanitize=False, removeHs=False):
    mols = list(Chem.SDMolSupplier(sdf_file, sanitize=sanitize, removeHs=removeHs))
    return mols[0] if len(mols) == 1 else mols

def write_sdf(sdf_file, mols):
    w = Chem.SDWriter(sdf_file)
    for mol in (mols if isinstance(mols, Iterable) else [mols]):
            w.write(mol)
    w.close()

def read_smi(smi_file, delimiter='\n', sanitize=False):
    mols = Chem.SmilesMolSupplier(smi_file, delimiter=delimiter, sanitize=sanitize)
    return mols

def to_mols(smiles):
    return [Chem.MolFromSmiles(smile) for smile in smiles]

def to_smiles(mols):
    return [Chem.MolToSmiles(mol) for mol in mols]

def standard_mol(mol):
    """Reset the molecule to *standard* form without conformation.  
    Not same as `Chem.MolStandardize.rdMolStandardize.Cleanup((Mol)`
    """
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))

def _sanitize_valid(mol):
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

@contextmanager
def temp_manager(suffix: str, output_dir: str):
    with tempfile.NamedTemporaryFile(
        suffix=suffix,
        dir=output_dir,
        delete=False
    ) as tmp:
        tmp_file = tmp.name
    try:
        yield tmp_file
    finally:
        try:
            os.unlink(tmp_file)
        except FileNotFoundError:
            pass


class Preprocess():
    def __init__(self, mols):
        self.mols = mols
        self.mols_num = len(mols)
        
    def valid(self):
        valids = [mol for mol in self.mols if _sanitize_valid(mol)]
        project_logger.info(f'Valid SMILES: {len(valids)} out of {self.mols_num}')
        return valids
    
    def unique(self, if_valid=True) -> Tuple[List[str], List[Chem.Mol]]:
        """Check the uniqueness of the molecule list.   
        It will return a deduplicated list of canonical SMILES `unique_smis` and their corresponding original Mol object lists `unique_mols`.  
        It is important to note that if there are different 3D conformations with the same canonical SMILES,
        the Mol objects will be assessed based on RMSD to determine if they represent the same conformation.   
        If no duplicates are found, all conformations will be retained.

        Args:
            if_valid (bool, optional): Using molecules list passed the valid check. Defaults to True.
        """
        mols = self.mols if not if_valid else self.valid()
        unique_di = defaultdict(list)
        for i, mol in enumerate(mols):
            smi = Chem.MolToSmiles(mol)
            unique_di[smi].append(i)

        consider_3D = False
        unique_smis, unique_mols = [], []
        for smi, indices in unique_di.items():
            unique_smis.append(smi)
            if len(indices) > 1:
                if mols[indices[0]].GetNumConformers(): # Consider conformer uniqueness
                    consider_3D = True
                    indices = _deduplicate_3D(mols, indices)
            unique_mols.extend([mols[idx] for idx in indices])
        
        project_logger.info(f'Unique SMILES: {len(unique_smis)} out of {self.mols_num}')
        if consider_3D:
            project_logger.info(f'3D Conformation checked...')
            project_logger.info(f'Unique Conformation: {len(unique_mols)} out of {self.mols_num}')
        return unique_smis, unique_mols