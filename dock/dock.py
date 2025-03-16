'''
Author: Rui Qin
Date: 2025-03-15 13:52:13
LastEditTime: 2025-03-16 18:15:42
Description: 
'''
import argparse
from utils.io import *
from utils.docking import LigPrep
from utils.logger import project_logger
from utils.constant import TARGETS
from dock.docking_vina import VinaDock
from dock.docking_gnina import GninaDock
from pathlib import Path

def breakpoint_check(result_pkl):
    if result_pkl.exists():
        finnal_idx = read_pkl(result_pkl)[-1]['index']
        project_logger.warning(f"Docking results already exists.  Docking will proceed starting from the {finnal_idx}th molecule.")
        return finnal_idx
    return 0


def setup_arguments(parser: argparse.ArgumentParser):
    parser.add_argument('-p', '--path', type=str, help='Path to the folder where generated molecules for testing will be stored.')
    parser.add_argument('-m', '--method', type=str, default='gnina', help='Docking method to use (`gnina` or `vina`). Default is gnina.')
    parser.add_argument('--verbose', type=int, default=0, help='Verbosity level of the docking process.')
    parser.add_argument('--seed', type=int, default=0, help='Random seed for docking.')
    parser.add_argument('--exhaust', type=int, default=8, help='Exhaustiveness of docking.')
    parser.add_argument('--poses', type=int, default=1, help='Number of poses to generate.')
    parser.add_argument('--mode', type=str, default='dock', help='Docking mode (`dock` or `score_only`). Default is dock.')
    parser.add_argument("--optimize", action="store_true", help="Optimized the original 3D conformation if available.")
    parser.add_argument("--reset", action="store_true", help="Reset the original 3D conformation if available.")
    return parser

def execute(args):
    work_dir = Path(args.path)
    for target in TARGETS:
        target_dir = work_dir / target
        if not target_dir.exists():
            project_logger.warning(f"Target folder {target_dir} not found.")
            continue

        # Breakpoint check
        Path(target_dir/'results').mkdir(parents=True, exist_ok=True)
        result_pkl = target_dir/'results/docking_results.pkl'
        if breakpoint_check(result_pkl):
        
        for ligand in target_dir.iterdir():

            #TODO Ligand preparation

            #TODO Docking