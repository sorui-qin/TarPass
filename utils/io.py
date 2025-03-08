'''
Author: Rui Qin
Date: 2025-03-08 15:38:31
LastEditTime: 2025-03-08 16:43:21
Description: 
'''
from rdkit import Chem
from collections.abc import Iterable

def write_sdf(sdf_file, mols):
    with Chem.SDWriter(sdf_file) as w:
        for mol in (mols if isinstance(mols, Iterable) else [mols]):
            w.write(mol)

def read_sdf(sdf_file, sanitize=0):
    mols = list(Chem.SDMolSupplier(sdf_file, sanitize=sanitize))
    return mols[0] if len(mols) == 1 else mols

def read_smi(smi_file):
    with open(smi_file, 'r') as f:
        smi_list = f.read.split('\n')
        mols = [Chem.MolFromSmiles(smi) for smi in smi_list]
    return mols