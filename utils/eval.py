'''
Author: Rui Qin
Date: 2025-06-13 19:02:11
LastEditTime: 2025-06-28 15:53:41
Description: 
'''
from pathlib import Path
from typing import Literal, Union
from utils.io import load_json
import pandas as pd

def find_dockpkl(results_dir:str|Path, mode:Literal['dock','origin','score_only']) -> Path:
    if (m:=mode) == 'origin':
        m = 'score_only'
    patterns = ['gnina', 'vina']
    for p in patterns:
        docked_pkl_path = Path(results_dir) / f'{p}-{m}_docking_results.pkl'
        if docked_pkl_path.exists():
            return docked_pkl_path
    raise FileNotFoundError(f"No docking results file found in '{results_dir}' for mode '{mode}'.\
                            Please rerun `tarpass dock -m {mode}` with a valid docking method.")

EVAL_MAP = {
    'dock': (load_json, 'dock_eval_results.json'),
    'mole': (pd.read_csv, 'mole_eval_results.csv')
}

def read_eval(results_dir:str|Path, mode:Literal['dock', 'mole']) -> Union[dict, pd.DataFrame]:
    """
    Read evaluation results from the specified directory based on the mode.
    """
    fn, fname = EVAL_MAP.get(mode, (None, None))
    if fn is None or fname is None:
        raise ValueError(f"Invalid mode '{mode}'. Expected one of {list(EVAL_MAP.keys())}.")
    eval_path = Path(results_dir) / fname
    if not eval_path.exists():
        raise FileNotFoundError(f"Evaluation file '{eval_path}' does not exist, please rerun `tarpass {mode}eval.")
    return fn(eval_path)