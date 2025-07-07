'''
Author: Rui Qin
Date: 2025-06-13 19:02:11
LastEditTime: 2025-07-06 16:19:16
Description: 
'''
from pathlib import Path
from typing import Literal, Union, overload
from utils.io import load_json
from utils.preprocess import to_inchi
import pandas as pd

def check_mols_equal(mols1, mols2):
    """Check if two sets of molecules are equal."""
    if len(mols1) != len(mols2):
        return False
    inchi1 = set(to_inchi(mols1))
    inchi2 = set(to_inchi(mols2))
    return inchi1 == inchi2

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

@overload
def read_eval(results_dir: str | Path, mode: Literal['dock']) -> list[dict]:
    ...

@overload
def read_eval(results_dir: str | Path, mode: Literal['mole']) -> pd.DataFrame:
    ...

def read_eval(results_dir:str|Path, mode:Literal['dock', 'mole']) -> Union[list[dict], pd.DataFrame]:
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