'''
Author: Rui Qin
Date: 2025-07-11 17:20:24
LastEditTime: 2025-07-11 20:34:02
Description: 
'''
import numpy as np
import pandas as pd
from analysis.collect_eval import AnalysisBase, DataCroupier
from module.similarity import Similarity

class MolDisAnalysis(AnalysisBase):
    """Analysis for molecular distance."""
    def __init__(self, collections: DataCroupier):
        super().__init__(collections, 'Mol')
        self.smis_croupier = [self._split_attr(attr) for attr in self.croupier_keys]
    # TODO: Lots of work to do...