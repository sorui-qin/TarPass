'''
Author: Rui Qin
Date: 2025-07-11 17:20:24
LastEditTime: 2025-07-17 18:10:53
Description: 
'''
import pandas as pd
from analysis.collect_eval import AnalysisBase, DataCroupier
from module.similarity import Similarity

class MolDisAnalysis(AnalysisBase):
    """Analysis for molecular distance."""
    def __init__(self, collections: DataCroupier):
        super().__init__(collections, 'Mol')
        self.smis_croupier = [self._split_attr(attr) for attr in self.croupier_keys]
        self.test_smis = self.smis_croupier[0]
    
    def _internal_distance(self) -> pd.DataFrame:
        simi = Similarity(self.test_smis) # type: ignore
        return pd.DataFrame([{
            'IntDiv': round(simi.intdiv(), 4),
            'IntDiv-Scaffold': round(simi.intdiv(cal_scaffold=True), 4),
            '#Circle': round(simi.circle(), 4),
        }])
        
    def _compare_distance(self, idx:int) -> pd.DataFrame:
        key_map = {
            0: 'Ref',
            1: 'Decoy'
        }
        key = key_map[idx]
        simi = Similarity(self.test_smis, self.smis_croupier[idx+1]) # type: ignore
        return pd.DataFrame([{
            f'{key}_similarity': round(simi.similarity(), 4),
            f'{key}_scaff_similarity': round(simi.scaffold_similarity(), 4),
            f'{key}_FCD': round(simi.fcd(), 4)
        }])
    
    def analysis(self):
        dfs = (
            [pd.DataFrame({'target': self.target}, index=[0])] +
            [self._internal_distance()] +
            [self._compare_distance(idx) for idx in range(2)]
        )
        return pd.concat(dfs, axis=1)

# TODO: Cross-targets analysis