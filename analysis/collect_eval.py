'''
Author: Rui Qin
Date: 2025-07-07 17:25:29
LastEditTime: 2025-07-07 20:45:13
Description: 
'''
from collections import namedtuple
from typing import Optional
from utils.constant import TARGETS, Path
from utils.io import read_pkl, write_pkl
from utils.eval import find_dockpkl, read_eval
from utils.preprocess import conformation_check
import pandas as pd

IS_3D = True

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
    
    def _get_mols(self) -> list:
        dockpkl = find_dockpkl(self.dir, mode='dock')
        return [res['mol'] for res in read_pkl(dockpkl)]
    
    def _read_eval(self) -> tuple[list[dict], pd.DataFrame]:
        dock_eval = read_eval(self.dir, mode='dock')
        mole_eval = read_eval(self.dir, mode='mole')
        if len(dock_eval) != len(mole_eval):
            raise ValueError(f"Docking and Molecule evaluation results length mismatch for {target}.")
        return dock_eval, mole_eval


for target in TARGETS:
    # TODO: 加上具体路径
    results_dir = Path() / f'{target}/results'
    # Read Mol
    dockpkl = find_dockpkl(results_dir, mode='dock')
    mols = [res['mol'] for res in read_pkl(dockpkl)]
    IS_3D = conformation_check(mols)

    dock_eval = read_eval(results_dir, mode='dock')
    mole_eval = read_eval(results_dir, mode='mole')
    if len(dock_eval) != len(mole_eval):
        raise ValueError(f"Docking and Molecule evaluation results length mismatch for {target}.")
    
    # Read Dock
    dock_dicts = [res['dock'] for res in dock_eval]
    score_dicts = [res['score'] for res in dock_eval] if IS_3D else []