'''
Author: Rui Qin
Date: 2025-07-10 15:07:11
LastEditTime: 2025-07-12 15:27:58
Description: 
'''
import numpy as np
import pandas as pd
from analysis.collect_eval import AnalysisBase, DataCroupier
from module.significance import SignificanceTester


def get_median_iqr(data):
    """Get median, IQR (Q1 & Q3) for the given data.
    Args:
        data (ArrayLike): Input data array.
    """
    data = np.array(data)[~np.isnan(data)]
    return (np.median(data),
            np.percentile(data, 25),
            np.percentile(data, 75))


class PLIAnalysis(AnalysisBase):
    """Analysis for protein-ligand interactions (PLI) related metrics from docking result."""
    def __init__(self, collections: DataCroupier):
        super().__init__(collections, 'Dock')
        
        self.numericals, self.interactions = [], []
        for attr in self.croupier_keys:
            num, inter = self._split_attr(attr) # type: ignore
            self.numericals.append(num)
            self.interactions.append(inter)
        self.test_numerical = self.numericals[0]
        self.test_interaction = self.interactions[0]
    
    def _get_median_iqr(self, key: str) -> pd.DataFrame:
        med, q1, q3 = get_median_iqr(self.test_numerical[key])
        return pd.DataFrame(
            {key: f'[{round(med, 4)}, {round(q1, 4)}, {round(q3, 4)}]'},
            index=[0])
    
    def _get_mean(self, key:str) -> pd.DataFrame:
        mean = np.nanmean(self.test_numerical[key])
        return pd.DataFrame({key: round(mean, 4)}, index=[0])
    
    def _score_significance(self, key='score') -> pd.DataFrame:
        control_data = self.test_numerical[key]
        data_groups = [i[key] for i in self.numericals[1:]]
        tester = SignificanceTester(data_groups, control_data, key)
        return tester.ref_decoy_analysis(alternative='less')
    
    def _count_interactions(self):
        total = len(self.test_interaction['detected_interactions'])
        all_inters = []
        for idx in range(total):
            interaction = self.test_interaction['detected_interactions'][idx]
            if interaction is None:
                all_inters.append(0)
            for v in interaction.values():
                all_inters.append(sum(len(i) for i in v))
        
        return pd.DataFrame({
            'average_interactions': round(np.mean(all_inters), 4),
        }, index=[0])


    def analysis(self) -> pd.DataFrame:
        affinity_cols = ['score', 'sucos', 'centroid shift', 'LE_heavyatom', 'LE_mw']
        boolean_cols = ['no_clashes', 'fully_matched', 'matched_rate']
        
        dfs = (
        [self._get_median_iqr(col) for col in affinity_cols] +
        [self._get_mean(col) for col in boolean_cols] +
        [self._score_significance()] +
        [self._count_interactions()]
    )
        return pd.concat(dfs, axis=1)


class PLIAnalysisScore(PLIAnalysis):
    """Analysis for protein-ligand interactions (PLI) related metrics from scoring result of 3D *in-situ* methods.
    """
    def __init__(self, collections: DataCroupier):
        super().__init__(collections)
        self.test = collections.test_data.Score
        self.test_numerical, self.test_interaction = self._split_attr('test') # type: ignore