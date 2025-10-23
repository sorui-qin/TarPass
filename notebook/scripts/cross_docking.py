'''
Author: Rui Qin
Date: 2025-10-23 14:17:35
LastEditTime: 2025-10-23 15:20:13
Description: Execute and anlysis cross-docking for two specific targets. 
It must be run after finishing the self-docking.
'''
import argparse
from pathlib import Path
from dock.dock_class import easydock
from utils.io import write_pkl, read_pkl, dump_json
from utils.logger import project_logger
from utils.constant import ROOT
from eval.dockeval import DockEval

def cross_dock_eval(b_dir, target_a):
    results_dir = b_dir / 'crossdocking'
    results_dir.mkdir(exist_ok=True)
    b_mols = [i['mol'] for i in read_pkl(b_dir / 'results/gnina-dock_docking_results.pkl')]
    results = easydock(b_mols, target_a)
    write_pkl(results_dir / f'crossdock_{target_a}.pkl', results)
    project_logger.info(f"Cross-docking to {target_a} completed.")
    
    poses = [res[0] for res in results]
    scores = [res[1] for res in results]
    eval_results = DockEval(poses, scores, target_a).evaluate()
    dump_json(results_dir / f'crossdock_{target_a}_eval.json', eval_results)
    project_logger.info(f"Cross-docking evaluation to {target_a} completed.")


def main(args):
    work_dir = Path(args.path)
    a_dir = work_dir / args.target_a
    b_dir = work_dir / args.target_b

    if not a_dir.exists() or not b_dir.exists():
        project_logger.error("One or both target directories do not exist.")
        return

    # A -> B
    cross_dock_eval(a_dir, args.target_b)
    # B -> A
    cross_dock_eval(b_dir, args.target_a)
    project_logger.info(f"Cross-docking completed.")

parser = argparse.ArgumentParser(description="Execute and analyze cross-docking for two specific targets.")
parser.add_argument('--path', '-p', type=str, required=True, help='Path to the docking results directory.')
parser.add_argument('--target_a', '-a', type=str, required=True, help='First target name for cross-docking.')
parser.add_argument('--target_b', '-b', type=str, required=True, help='Second target name for cross-docking.')
args = parser.parse_args()
main(args)