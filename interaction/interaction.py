'''
Author: Rui Qin
Date: 2025-05-13 23:35:54
LastEditTime: 2025-06-13 17:11:04
Description: 
'''
import argparse
import json
from typing import Literal
from interaction.interaction_tools import Path, interactions, plip_tmp
from utils.constant import DASHLINE, ROOT, TARGETS
from utils.io import append_pkl, read_pkl
from utils.logger import log_config, project_logger
from utils.preprocess import read_in


def find_docked_pkl(results_dir:Path, mode:Literal['dock','origin']) -> Path|None:
    if (m:=mode) == 'origin':
        m = 'score_only'
    patterns = ['gnina', 'vina']
    for p in patterns:
        docked_pkl_path = results_dir / f'{p}-{m}_docking_results.pkl'
        if docked_pkl_path.exists():
            return docked_pkl_path
    return None

def setup_arguments(parser: argparse.ArgumentParser):
    parser.add_argument('--mode', required=True, type=str, choices=['dock', 'origin'], help='Anlysis interactions with docking pose or original generated pose.')
    return parser

def execute(args):
    log_config(project_logger, args)
    # Preparation for reading poses
    work_dir = Path(args.path)
    mode = args.mode
    # Analyze
    allkey_inters = json.load(open(ROOT/'interaction/key_interactions.json'))

    for target in TARGETS:
        plip_tmp() # Clean up temp dir
        project_logger.info(f"Analyzing interactions with {target}...")
        key_inters = allkey_inters[target]
        empty_pdb = next((ROOT / f'Targets/{target}').glob(f'{target}*.pdb'))
        
        # Read in the generated molecules
        target_dir = work_dir / target
        if not target_dir.exists():
            project_logger.error(f"Target directory {target_dir} does not exist.")
            continue
        results_dir = target_dir / 'results'
        results_dir.mkdir(parents=True, exist_ok=True)
        _, mols = read_in(target_dir)
        
        if mode == 'dock':
            docked_pkl = find_docked_pkl(results_dir, mode=mode)
            if not docked_pkl:
                raise ValueError(f"Docking results not found in {target_dir}, please running `tarpass dock`.")
            docked_results = read_pkl(docked_pkl)
            if len(docked_results) < len(mols):
                raise ValueError(f"Docking results are incomplete, please rerun `tarpass dock`.")
            # Get the docking poses
            poses = [docked['pose'] for docked in docked_results]
        elif mode == 'origin':
            poses = mols

        # Check if have already analyzed
        store_pkl = results_dir / f'{mode}_interactions.pkl'
        if store_pkl.exists():
            processed = read_pkl(store_pkl)
            if len(processed) == len(poses):
                project_logger.info(f"Interaction analysis already completed for {len(poses)} molecules.")
                continue
        
        # Analyze interactions
        all_match = interactions(poses, empty_pdb, key_inters)

        # Save results
        store_li = []
        for mol, match in zip(mols, all_match):
            store_li.append({**{'idx': int(mol.GetProp('_Name'))}, **match})
        append_pkl(store_pkl, store_li)
        project_logger.info(f"Interaction analysis result is save at {results_dir}.")
        project_logger.info(DASHLINE)