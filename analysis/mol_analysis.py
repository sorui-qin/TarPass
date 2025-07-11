'''
Author: Rui Qin
Date: 2025-07-11 17:20:24
LastEditTime: 2025-07-11 17:35:54
Description: 
'''
import numpy as np
import pandas as pd
from analysis.collect_eval import AnalysisBase, DataCroupier
import fcd

class MolAnalysis(AnalysisBase):
    """Analysis for molecular properties and drug-likeness metrics."""
    def __init__(self, collections: DataCroupier):
        super().__init__(collections, 'Mol')