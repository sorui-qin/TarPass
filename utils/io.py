'''
Author: Rui Qin
Date: 2025-03-08 15:38:31
LastEditTime: 2025-03-17 16:44:14
Description: 
'''
import os
import pickle
import tempfile
from rdkit import Chem
from collections.abc import Iterable
from contextlib import contextmanager
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def write_pkl(pkl_file, data):
    with open(pkl_file, 'wb') as f:
        pickle.dump(data, f)

def append_pkl(pkl_file, data):
    with open(pkl_file, 'ab') as f:
        pickle.dump(data, f)

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

def read_smi(smi_file, delimiter='', sanitize=False):
    supp = Chem.SmilesMolSupplier(smi_file, delimiter=delimiter, sanitize=sanitize, titleLine=False)
    mols = [mol for mol in supp if mol]
    return mols[0] if len(mols) == 1 else mols

@contextmanager
def temp_manager(suffix: str, output_dir='./tmp', auto_remove=True):
    os.makedirs(output_dir, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        suffix=suffix,
        dir=output_dir,
        delete=False
    ) as tmp:
        tmp_file = tmp.name
    try:
        yield tmp_file
    finally:
        if auto_remove:
            try:
                os.unlink(tmp_file)
            except FileNotFoundError:
                pass