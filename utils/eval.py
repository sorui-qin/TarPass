'''
Author: Rui Qin
Date: 2025-06-13 19:02:11
LastEditTime: 2025-06-13 19:02:14
Description: 
'''
from pathlib import Path
from typing import Literal

def find_dockpkl(results_dir:Path, mode:Literal['dock','origin','score_only']) -> Path|None:
    if (m:=mode) == 'origin':
        m = 'score_only'
    patterns = ['gnina', 'vina']
    for p in patterns:
        docked_pkl_path = results_dir / f'{p}-{m}_docking_results.pkl'
        if docked_pkl_path.exists():
            return docked_pkl_path
    return None