'''
Author: Rui Qin
Date: 2025-07-07 17:25:29
LastEditTime: 2025-07-08 15:55:35
Description: 
'''
from collections import namedtuple
from typing import Optional, Literal
from utils.constant import TARGETS, Path
from utils.io import read_pkl, write_pkl
from utils.eval import find_dockpkl, read_eval
from utils.preprocess import conformation_check
import pandas as pd


class CollectEval:
    """
    Collect evaluation results for docking and molecule evaluation.
    """
    def __init__(self, results_dir: Path):
        self.dir = results_dir
        self.target = self.dir.parent.name
        self.mols = self._get_mols()
        self.confs = conformation_check(self.mols)
        self.dock_eval, self.mole_eval = self._read_eval()
        self.dock_dicts = [res['dock'] for res in self.dock_eval]
        self.score_dicts = [res['score'] for res in self.dock_eval] if self.confs else []
    
    def _get_mols(self) -> list:
        dockpkl = find_dockpkl(self.dir, mode='dock')
        return [res['mol'] for res in read_pkl(dockpkl)]
    
    def _read_eval(self) -> tuple[list[dict], pd.DataFrame]:
        dock_eval = read_eval(self.dir, mode='dock')
        mole_eval = read_eval(self.dir, mode='mole')
        if len(dock_eval) != len(mole_eval):
            raise ValueError(f"Docking and Molecule evaluation results length mismatch for {self.target}.")
        return dock_eval, mole_eval
    
    def _split_dicts(self, mode=Literal['dock', 'score'], split_key:str='matched_rate'):
        dicts = self.dock_dicts if mode == 'dock' else self.score_dicts
        if not dicts:
            raise ValueError(f"No evaluation results found for mode '{mode}' in target '{self.target}'.")
        
        numericals, interactions = [], []
        for d in dicts:
            if split_key not in d:
                raise ValueError(f"Split key '{split_key}' not found in dictionary.")
            
            items = list(d.items())

            # Default split index for matched_rate
            split_idx = None if split_key != 'matched_rate' else 10
            for i, (key, _) in enumerate(items):
                if key == split_key:
                    split_idx = i + 1
                    break
            
            numericals.append(dict(items[:split_idx]))
            interactions.append(dict(items[split_idx:]))
        return numericals, interactions