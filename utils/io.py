'''
Author: Rui Qin
Date: 2025-03-08 15:38:31
LastEditTime: 2025-09-10 11:59:11
Description: 
'''
import json
import os
import pickle
import shutil
import tempfile
from collections.abc import Iterable
from contextlib import contextmanager
from typing import Any
import yaml
from rdkit import Chem, RDLogger
from utils.constant import Path, TMP
from utils.logger import project_logger

RDLogger.DisableLog('rdApp.*') # type: ignore

def write_pkl(pkl_file, data):
    with open(pkl_file, 'wb') as f:
        pickle.dump(data, f)

def append_pkl(pkl_file, data):
    with open(pkl_file, 'ab') as f:
        pickle.dump(data, f)

def read_pkl(pkl_file) -> Any:
    data = []
    with open(pkl_file, 'rb') as f:
        while True:
            try:
                obj = pickle.load(f)
                data.append(obj)
            except EOFError:
                break
    return data if len(data) > 1 else data[0]

def load_json(json_file):
    return json.load(open(json_file, 'r'))

def dump_json(json_file, data):
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)

def read_pdb_rdmol(pdb_file:str|Path, sanitize=False, removeHs=False):
    mol = Chem.MolFromPDBFile(str(pdb_file), 
                              sanitize=sanitize, removeHs=removeHs)
    return mol

def read_sdf(sdf_file:str|Path, sanitize=False, removeHs=False):
    mols = [mol for mol in Chem.SDMolSupplier(
        str(sdf_file), sanitize=sanitize, removeHs=removeHs) if mol]
    return mols[0] if len(mols) == 1 else mols

def write_sdf(sdf_file, mols):
    w = Chem.SDWriter(sdf_file)
    for mol in (mols if isinstance(mols, Iterable) else [mols]):
        w.write(mol)
    w.close()

def read_smi(smi_file:str|Path, delimiter='', sanitize=False):
    supp = Chem.SmilesMolSupplier(str(smi_file), 
                                  delimiter=delimiter, sanitize=sanitize, titleLine=False)
    mols = [mol for mol in supp if mol]
    return mols[0] if len(mols) == 1 else mols

def read_strings(string_file):
    if str(string_file).split('.')[-1] not in ['smi', 'txt']:
        project_logger.warning(f"File {string_file} is not a valid string file.")
        return []
    with open(string_file, 'r') as f:
        lines = [line.strip() for line in f]
    return lines

def read_yaml(yaml_file):
    with open(yaml_file, 'r') as f:
        data = yaml.safe_load(f)
    return data

def temp_dir(output_dir=TMP):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

@contextmanager
def temp_manager(suffix:str, output_dir=TMP, auto_remove=True):
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