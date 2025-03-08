'''
Author: Rui Qin
Date: 2024-01-04 17:22:05
LastEditTime: 2025-03-08 20:47:35
Description: 
'''
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
import numpy as np
from chemmeasure.measures import NCircles
from utils.measure import *


def define_measures(fps):
    vectorizer = fingerprints
    sim_mat_func = similarity_matrix_tanimoto
    measure = NCircles(vectorizer=vectorizer, sim_mat_func=sim_mat_func, threshold=0.75)
    circ, _ = measure.measure(fps, is_vec=True, n_chunk=64)
    return circ


def similarity(test_mols, standard_mols):
    """
    Calculate the similarity between two molecular sets.
    Input:
        test_mols : Set to be tested, typically the generated set.
        standard_mols : Set to be compared, such as a so-called 'test set' or ground truth set.
    """
    matrix = similarity_matrix_tanimoto(
        fingerprints(test_mols), 
        fingerprints(standard_mols))
    return matrix.sum() / (len(test_mols) * len(standard_mols))

if __name__ == '__main__':
    from utils.io import read_sdf
    mols = read_sdf('testfile/5HT2A_gen.sdf')
    fps = fingerprints(mols)
    print(define_measures(fps))