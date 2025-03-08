'''
Author: Rui Qin
Date: 2024-01-04 17:22:05
LastEditTime: 2025-03-09 01:19:44
Description: 
'''
from rdkit.Chem.Scaffolds import MurckoScaffold
import numpy as np
from chemmeasure.measures import NCircles
from utils.measure import *

class Similarity():
    def __init__(self, mols:list, ref_mols=None, compared_mode=False):
        self.mols = mols
        self.fps = fingerprints(self.mols)
        self.compared = compared_mode
        if compared_mode:
            if not ref_mols:
                raise ValueError('Under comparison mode, please provide reference molecules set')
            self.ref = ref_mols


    def _not_compared(self):
        if self.compared:
            raise ValueError('This function is not suitable under comparison mode!')


    def circle(self, thershold=0.75) -> float:
        """Calculate #Circle. 

        Args:
            thershold (float, optional): Thershold of #Circle. Defaults to 0.75.

        Returns:
            float: Vaule of #Circle.
        """
        self._not_compared()
        vectorizer = fingerprints
        sim_mat_func = similarity_matrix_tanimoto
        measure = NCircles(vectorizer=vectorizer, sim_mat_func=sim_mat_func, threshold=thershold)
        circ, _ = measure.measure(self.fps, is_vec=True, n_chunk=64)
        return circ


    def similarity(self, limit=10000, seed=0) -> float:
        """Calculate Tanimoto similarity with ECFP4.

        Args:
            limit (int, optional): Maximum number of molecules for calculation. \n
            If exceeded, molecules will be randomly selected up to the limit. Defaults to 10000.
            seed (int, optional): Random seed. Defaults to 0 (not set).

        Returns:
            float: Value of Tanimoto similarity
        """
        fps1 = sample(self.fps, n=limit, fixed_seed=seed)
        fps2 = fps1 if self.compared else fingerprints(sample(self.ref, n=limit, fixed_seed=seed))
        matrix = similarity_matrix_tanimoto(fps1, fps2)
        return matrix.sum() / (len(fps1) * len(fps2))