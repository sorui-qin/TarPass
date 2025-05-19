'''
Author: Rui Qin
Date: 2025-05-13 23:35:54
LastEditTime: 2025-05-19 20:23:55
Description: 
'''
import argparse
import json
from tqdm import tqdm
from interaction.interaction_tools import *
from multiprocessing import Pool
from functools import partial
from utils.io import read_pkl, append_pkl
from utils.preprocess import read_in
from utils.logger import project_logger, log_config
from utils.constant import ROOT, TARGETS, DASHLINE


def find_docked_pkl(results_dir:Path) -> Path|None:
    patterns = ['gnina', 'vina']
    for p in patterns:
        docked_pkl_path = results_dir / f'{p}-dock_docking_results.pkl'
        if docked_pkl_path.exists():
            return docked_pkl_path
    return None

def setup_arguments(parser: argparse.ArgumentParser):
    #parser.add_argument('-p', '--path', required=True, type=str, help='path to the folder where generated molecules for testing will be stored.')
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
            docked_pkl = find_docked_pkl(results_dir)
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
            else:
                project_logger.info(f"Detected previous interaction analysis results, resuming from {len(processed)} of {len(poses)} molecules.")
                poses = poses[len(processed):]

        # Check if all molecules are 3D conformations
        if not all([pose.GetNumConformers() for pose in poses]):
            raise RuntimeError(f"Some molecules have not 3D conformations. Please check the input molecules.")
        total_num = len(poses)
        project_logger.info(f"Total {total_num} molecules to analyze.")

        # Running interaction analysis in parallel
        analyze_func = partial(analyze_tmppdb, empty_pdb=empty_pdb, key_inters=key_inters)
        with Pool() as pool:
            all_match = list(tqdm(pool.imap(analyze_func, poses),
                                desc="Analyzing interactions", total=total_num))
        # Save results
        store_li = []
        for mol, match in zip(mols, all_match):
            store_li.append({**{'idx': mol.GetProp('_Name')}, **match})
        append_pkl(store_pkl, store_li)
        project_logger.info(f"Interaction analysis result is save at {results_dir}.")
        project_logger.info(DASHLINE)