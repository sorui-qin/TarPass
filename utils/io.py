'''
Author: Rui Qin
Date: 2025-03-08 15:38:31
LastEditTime: 2025-03-10 11:05:51
Description: 
'''
from rdkit import Chem
from collections.abc import Iterable

def write_sdf(sdf_file, mols):
    with Chem.SDWriter(sdf_file) as w:
        for mol in (mols if isinstance(mols, Iterable) else [mols]):
            w.write(mol)

def read_sdf(sdf_file, sanitize=False, removeHs=False):
    mols = list(Chem.SDMolSupplier(sdf_file, sanitize=sanitize, removeHs=removeHs))
    return mols[0] if len(mols) == 1 else mols

def read_smi(smi_file):
    with open(smi_file, 'r') as f:
        smi_list = f.read.split('\n')
        mols = [Chem.MolFromSmiles(smi) for smi in smi_list]
    return mols

def standard_mol(mol):
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))