'''
Author: Rui Qin
Date: 2025-03-08 15:00:12
LastEditTime: 2025-06-18 19:02:10
Description: 
'''
# Partly adapted from https://github.com/yutxie/chem-measure/blob/main/utils.py
import random
import numpy as np
from rdkit.Chem import DataStructs, Mol
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from multiprocessing import Pool


def sample(anylist, n=10000, fixed_seed=0):
    if fixed_seed: random.seed(42)
    if len(anylist) > n:
        return random.sample(anylist, n)
    else:
        return anylist
    
mfpgen = GetMorganGenerator(radius=2, fpSize=2048) # ECFP4

fp_types = {
    'fp': mfpgen.GetFingerprint,
    'fpnp': mfpgen.GetFingerprintAsNumPy,
    'sfp': mfpgen.GetSparseFingerprint,
    'cfp': mfpgen.GetCountFingerprint,
    'scfp': mfpgen.GetSparseCountFingerprint
}

def fingerprint(mol:Mol, fp_type:str):
    return fp_types[fp_type](mol)

def fingerprints(mols:list[Mol], fp_type='fp'):
    return [fingerprint(mol, fp_type) for mol in mols]

def morgan_frags(mol:Mol):
    fp = fingerprint(mol, fp_type='scfp')
    return list(fp.GetNonzeroElements().keys())

def similarity(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def similarities_tanimoto(fp, fps):
    return DataStructs.BulkTanimotoSimilarity(fp, fps)

def similarity_matrix_tanimoto(fps1, fps2):
    similarities = [DataStructs.BulkTanimotoSimilarity(fp, fps2) for fp in fps1]
    return np.array(similarities)

def similarity_matrix_parallel(fps1, fps2, chunksize=100):
    from functools import partial
    similar_func = partial(similarities_tanimoto, fps=fps2)
    
    pool = Pool()
    similarities = list(pool.imap(similar_func, fps1, chunksize=chunksize))
    return np.array(similarities)