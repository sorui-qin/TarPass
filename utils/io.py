'''
Author: Rui Qin
Date: 2025-03-08 15:38:31
LastEditTime: 2025-03-13 16:50:49
Description: 
'''
import pickle
from rdkit import Chem
from collections.abc import Iterable
from utils.logger import project_logger

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
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))


class Preprocess():
    def __init__(self, mols):
        self.mols = mols
        self.mols_num = len(mols)
        self.smiles = to_smiles(mols)
        
    def valid(self):
        valids = list(filter(None, to_mols(self.smiles)))
        project_logger.info(f'Valid SMILES: {len(valids)} out of {self.mols_num}')
        return valids
    
    def unique(self, if_valid=True):
        smiles = self.smiles if not if_valid else to_smiles(self.valid())
        smiles = list(set(smiles))
        project_logger.info(f'Unique SMILES: {len(smiles)} out of {self.mols_num}')
        return to_mols(smiles)