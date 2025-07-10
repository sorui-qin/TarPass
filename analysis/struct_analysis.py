'''
Author: Rui Qin
Date: 2025-07-10 15:07:11
LastEditTime: 2025-07-10 20:55:47
Description: 
'''
from typing import Literal, Optional
from analysis.collect_eval import MoleculesData, StructInfo, AnalysisBase, DataCroupier
import numpy as np
import pandas as pd
from typing import cast
from module.significance import SignificanceTester


def get_median_iqr(data):
    """Get median, IQR (Q1 & Q3) for the given data.
    Args:
        data (ArrayLike): Input data array.
    """
    return (np.median(data),
            np.percentile(data, 25),
            np.percentile(data, 75))


class StructAnalysis(AnalysisBase):
    def __init__(self, collections: DataCroupier):
        super().__init__(collections, 'Dock')
        
        self.numericals, self.interactions = [], []
        for attr in self.croupier:
            num, inter = self._split_attr(attr) # type: ignore
            self.numericals.append(num)
            self.interactions.append(inter)
        self.test_numerical = self.numericals[0]
        self.test_interaction = self.interactions[0]
    
    def _get_median_iqr(self, key: str):
        med, q1, q3 = get_median_iqr(self.test_numerical[key])
        return pd.DataFrame(
            {key: f'[{round(med, 4)}, {round(q1, 4)}, {round(q3, 4)}]'},
            index=[0]
            )
    
    def _score_significance(self, key='score'):
        control_data = self.test_numerical[key]
        data_groups = [i[key] for i in self.numericals[1:]]
        tester = SignificanceTester(data_groups, control_data, key)
        return tester.ref_decoy_analysis(alternative='less')