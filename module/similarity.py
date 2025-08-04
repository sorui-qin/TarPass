'''
Author: Rui Qin
Date: 2024-01-04 17:22:05
LastEditTime: 2025-08-04 21:12:03
Description: 
'''
from typing import Optional
from fcd import get_fcd
from rdkit.Chem import Mol
from rdkit.Chem.Scaffolds import MurckoScaffold  # type: ignore[import-untyped]
from module.chemmeasure.measures import NCircles
from utils.measure import (fingerprints, sample, similarity_matrix_parallel,
                           similarity_matrix_tanimoto)
from utils.preprocess import to_mols


def scaffolds(mols):
    return [MurckoScaffold.GetScaffoldForMol(mol) for mol in mols]


def get_similarity(mols1:list[Mol], mols2:Optional[list[Mol]], 
                   limit=10000, seed=0, use_ref:bool=True) -> float:
    """Calculate Tanimoto similarity between two sets of molecules.
    Args:
        mols1 (list[Mol]): List of RDKit molecules for the first set.
        mols2 (Optional[list[Mol]]): List of RDKit molecules for the second set.
        limit (int, optional): Maximum number of molecules for calculation.  
            If exceeded, molecules will be randomly selected up to the limit. Defaults to 10000.
        seed (int, optional): Random seed for reproducibility. Defaults to 0 (not set).
        use_ref (bool, optional): Whether to use the second set of molecules as a reference. 
            If False, the first set will be used for both calculations. Defaults to True.
    """
    mols1 = sample(mols1, n=limit, fixed_seed=seed)
    fps1 = fingerprints(mols1, fp_type='fp') # Use bit vector for Tanimoto similarity
    if mols2 and use_ref:
        mols2 = sample(mols2, n=limit, fixed_seed=seed)
        fps2 = fingerprints(mols2, fp_type='fp')
    else:
        fps2 = fps1
    # matrix = similarity_matrix_parallel(fps1, fps2)
    matrix = similarity_matrix_tanimoto(fps1, fps2)
    return matrix.sum() / (len(fps1) * len(fps2))


class Similarity:
    def __init__(self, smis:list[str], ref_smis:Optional[list[str]]=None, limit=10000, seed=0):
        self.smis = smis
        self.mols = to_mols(smis)
        self.ref_smis = ref_smis
        self.ref = to_mols(ref_smis) if ref_smis else None
        self.limit = limit
        self.seed = seed

    def circle(self, thershold=0.75) -> float:
        """Calculate #Circle. Please see Xie et al., *How Much Space Has Been Explored? 
        Measuring the Chemical Space Covered by Databases and Machine-Generated Molecules.* ICLR 2023.

        Args:
            thershold (float, optional): Thershold of #Circle. Defaults to 0.75.

        Returns:
            float: Vaule of #Circle.
        """
        vectorizer = fingerprints
        sim_mat_func = similarity_matrix_tanimoto
        measure = NCircles(vectorizer=vectorizer, sim_mat_func=sim_mat_func, threshold=thershold)
        fps = fingerprints(self.mols, fp_type='fp')
        circ, _ = measure.measure(fps, is_vec=True, n_chunk=64)
        return circ

    def similarity(self, use_ref:bool=True) -> float:
        return get_similarity(self.mols,
                              self.ref if use_ref else None,
                              limit=self.limit, seed=self.seed,
                              use_ref=use_ref)

    def scaffold_similarity(self, use_ref:bool=True) -> float:
        """Calculate Tanimoto similarity of scaffold with ECFP4.
        """
        scaff1 = scaffolds(sample(self.mols, n=self.limit, fixed_seed=self.seed))
        scaff2 = scaffolds(sample(self.ref, n=self.limit, fixed_seed=self.seed)) if self.ref else None
        return get_similarity(scaff1, scaff2, use_ref=use_ref)

    def intdiv(self, cal_scaffold=False):
        func = self.similarity if not cal_scaffold else self.scaffold_similarity
        return 1 - func(use_ref=False)

    def fcd(self):
        if not self.ref_smis:
            raise ValueError("Reference SMILES list is required for FCD calculation.")
        return get_fcd(self.smis, self.ref_smis)