'''
Author: Rui Qin
Date: 2025-05-13 23:35:54
LastEditTime: 2025-05-15 11:30:52
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
        results_dir = target_dir / 'results'
        results_dir.mkdir(parents=True, exist_ok=True)
        docked_pkl = find_docked_pkl(results_dir)
        if docked_pkl is None:
            if mode == 'dock':
                raise ValueError(f"Docking results not found in {target_dir}, please running `tarpass dock`.")
            elif mode == 'origin':
                project_logger.warning(f"Docking results not found in {target_dir}, searching for original poses.")
                _, poses = read_in(target_dir)
        else:
            docked_list = read_pkl(docked_pkl)
            sele = 'pose' if mode == 'dock' else 'mol'
            poses = [pose[sele] for pose in docked_list]

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
        for pose, match in zip(poses, all_match):
            store_li.append({**{'idx': pose.GetProp('_Name')}, **match})
        append_pkl(results_dir / f'interactions.pkl', store_li)
        project_logger.info(f"Interaction analysis result is save at {results_dir}.")
        project_logger.info(DASHLINE)