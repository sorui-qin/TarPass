'''
Author: Rui Qin
Date: 2024-01-04 17:22:05
LastEditTime: 2024-12-26 19:11:53
Description: 
'''
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.spatial.distance import cosine as cos_distance
#from multiprocessing import Pool
#from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
import random
import numpy as np
import re

def read_sdf(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file, sanitize=False)
    mols = [i for i in supp]
    return mols

def read_smi(smi_file):
    with open(smi_file, 'r') as f:
        smi_list = f.read.split('\n')
        mols = [Chem.MolFromSmiles(smi) for smi in smi_list]
    return mols

def sample(mols, n=10000):
    if len(mols) > n:
        random.seed(42)
        sample_batch = random.sample(mols, n)
    else:
        sample_batch = mols
    return sample_batch

def get_filename(filename):
    return re.search(r'[^\\/]+(?=\.)', filename).group(0)

def cosine_sim(vec1, vec2):
    cos_sim = 1 - cos_distance(vec1, vec2)
    return cos_sim

def sim(args):
    fp1, fp2 = args
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def calc_smi_parallel(test_fps, standard_fps):
    test = np.array(test_fps)
    std = np.array(standard_fps)
    pairs = np.array(np.meshgrid(test, std)).T.reshape(-1, 2)
    results = process_map(
            sim,
            pairs,
            chunksize=1000
        )
    return np.array(results, dtype=float)


def getsimilarity(test_mols, standard_mols):
    """
    Calculate the similarit
    y between two molecular sets.
    Input:
        test_mols : Set to be tested, typically the generated set.
        standard_mols : Set to be compared, such as a so-called 'test set' or ground truth set.
        But actually similarity is symmetric,
        so how you input the both molecular lists doesn't affect the final result XD
    """
    
    test_batch, standard_batch = sample(test_mols, 5000), sample(standard_mols)
    test_fps = [AllChem.GetMorganFingerprint(mol, 2) for mol in test_batch]
    standard_fps = [AllChem.GetMorganFingerprint(mol, 2) for mol in standard_batch]
    similar_array = calc_smi_parallel(test_fps, standard_fps)
    similarity = similar_array.sum() / (len(test_batch) * len(standard_batch))
    return similarity


def scaffold_similarity(test_mols, standard_mols):
    """
    Calculate the scaffold similarity between two molecular sets.
    Input:
        test_mols : Set to be tested, typically the generated set.
        standard_mols : Set to be compared, such as a so-called 'test set' or ground truth set.
        But actually similarity is symmetric,
        so how you input the both molecular lists doesn't affect the final result XD
    """
    
    test_batch, standard_batch = sample(test_mols), sample(standard_mols)
    test_fps = [AllChem.GetMorganFingerprint(MurckoScaffold.GetScaffoldForMol(mol), 2) for mol in test_batch]
    standard_fps = [AllChem.GetMorganFingerprint(MurckoScaffold.GetScaffoldForMol(mol), 2) for mol in standard_batch]
    similar_array = calc_smi_parallel(test_fps, standard_fps)
    similarity = similar_array.sum() / (len(test_batch) * len(standard_batch))
    return similarity



if __name__ == '__main__':
    """
    The input files here are flexible,
    please pay attention and check the format when you input,
    then change the corresponding codes.
    """
    test_file = '' # Path of the test file, sdf file suggested.
    standard_flie = '' # Path of the standard_file, smi/txt file with SMILES strings suggested.
    test_mols = read_sdf(test_file)
    standard_mols = read_smi(standard_flie)

    similarity = "%.3f" % getsimilarity(test_mols, standard_mols)
    print (f'Similarity between {get_filename(test_file)} and {get_filename(standard_flie)} is {similarity}')


