'''
Author: Rui Qin
Date: 2025-03-15 13:52:13
LastEditTime: 2025-08-08 16:36:58
Description: 
'''
import argparse
from tqdm import tqdm
from utils.io import read_pkl, append_pkl
from utils.logger import project_logger, log_config
from utils.constant import ROOT, TARGETS, DASHLINE
from utils.preprocess import conformation_check, read_in
from pathlib import Path
from dock.dock_class import SingleDock, BatchDock

def breakpoint_check(result_pkl: Path, total_lens:int) -> int:
    if result_pkl.exists():
        project_logger.info(f"Detected previous docking results.")
        results = read_pkl(result_pkl)
        if isinstance(results, list):
            latest_idx = len(results)
            if latest_idx > total_lens:
                raise ValueError(f"Index {latest_idx} exceeds the total number of ligands {total_lens}.")
            project_logger.info(f"Docking will start from {latest_idx+1} of {total_lens} molecules.")
            return latest_idx
    return 0

def save_results(pkl_path, index, mol, pose, score, mode):
    """Save docking results to a pickle file."""
    if mode == 'dock':
        append_pkl(pkl_path, {'index': index, 'mol': mol, 'pose': pose, 'docking score': score})
    elif mode == 'score_only':
        append_pkl(pkl_path, {'index': index, 'mol': mol, 'vina score': score})

############## Execution Functions ##############

def setup_arguments(parser: argparse.ArgumentParser):
    group1 = parser.add_argument_group("Docking parameters\ndefaultly loaded from config file `configs/dock/gnina_dock.yml`")
    group1.add_argument('-m', '--mode', type=str, choices=['dock', 'score_only', 'all'], default='all', help='docking mode.')
    group1.add_argument('--verbose', type=int, help='verbosity level of the docking process.')
    group1.add_argument('--seed', type=int, help='random seed for docking.')
    group1.add_argument('--exhaust', type=int, help='exhaustiveness of docking.')
    group1.add_argument('--poses', type=int, help='number of poses to generate.')
    group1.add_argument('--config', type=str, default=f'{ROOT}/configs/dock/gnina_dock.yml', help='path to the configuration file.')
    group1.add_argument('--method', type=str, choices=['gnina', 'vina'], help='docking method to use (default: gnina).')

    group2 = parser.add_argument_group("Optional arguments")
    group2.add_argument('--reset', action="store_true", help="reset the original 3D conformation if available.")
    group2.add_argument('--in_single', action="store_true", default=False, help="docking in single mode.")
    return parser

def execute(args):
    log_config(project_logger, args)
    work_dir = Path(args.path)
    for target in TARGETS:
        target_dir = work_dir / target
        if not target_dir.exists():
            project_logger.warning(f"Target folder {target_dir} not found.")
            continue
        if not list(target_dir.glob('*')):
            project_logger.warning(f"No files found in {target_dir}.")
            continue
        project_logger.info(f'Start performing docking in {target}...')

        # Preprocess
        _, mols = read_in(target_dir, args.num)
        total_lens = len(mols)

        if args.mode == 'all':
            if conformation_check(mols):
                modes = ['dock', 'score_only'] 
            else:
                project_logger.warning(f"Conformations are not detected in {target}, switching to `dock` mode.")
                modes = ['dock']
        else:
            modes = [args.mode]

        for mode in modes:
            # Set docking mode
            project_logger.info(f"Current docking mode: {mode}.")
            dock_args = argparse.Namespace(**vars(args))
            dock_args.mode = mode

            # Breakpoint check
            Path(target_dir/'results').mkdir(parents=True, exist_ok=True)
            if args.reset:
                result_pkl = target_dir/f'results/{args.method}-reset_docking_results.pkl'
            else:
                result_pkl = target_dir/f'results/{args.method}-{mode}_docking_results.pkl'
            latest_idx = breakpoint_check(result_pkl, total_lens)
            if latest_idx == total_lens:
                project_logger.info(f"All molecules in {target} have been docked.")
                continue
            
            # Docking
            if args.in_single or latest_idx != 0:
                for i, mol in tqdm(enumerate(latest:=mols[latest_idx:]), 
                                desc=f'Docking with {target}', total=len(latest)):
                    index = range(total_lens)[latest_idx+i]
                    dock = SingleDock(mol, target, dock_args)
                    pose, score = dock.run()
                    # Save results
                    save_results(result_pkl, index, mol, pose, score, mode)
            
            elif not (args.in_single and latest_idx):
                # Execute batch docking if not in single mode and no previous results
                project_logger.info(f'Batch docking with {target}, total {total_lens}...')
                dock = BatchDock(mols, target, dock_args)
                results = dock.run()
                # Save results
                for index, (pose, score) in enumerate(results):
                    save_results(result_pkl, index, mols[index], pose, score, mode)
                if len(results) != total_lens:
                    raise ValueError(f"Batch docking results length {len(results)} does not match total ligands {total_lens}.")
                
            project_logger.info(f'Docking in {target} completed. Results saved in {result_pkl}.')
            project_logger.info(DASHLINE)