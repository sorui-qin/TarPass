'''
Author: Rui Qin
Date: 2025-06-13 11:22:44
LastEditTime: 2025-06-14 17:20:27
Description: 
'''
from functools import partial
from pathlib import Path
from typing import Literal
import pandas as pd
from tqdm.contrib.concurrent import process_map
from eval.module.druglikeness import calc_le
from eval.module.sucos import check_sucos
from interaction.interaction_tools import allkey_inters, interactions
from utils.constant import DASHLINE, ROOT, TARGETS
from utils.eval import find_dockpkl
from utils.io import read_pkl, read_sdf
from utils.logger import log_config, project_logger
from utils.preprocess import Chem, conformation_check, read_in, to_smiles

TARGET_PATH = ROOT / 'Targets'

class DockEval():
    """Evaluate docking results of generated molecules for a specific target.
    Args:
        poses (list[Chem.Mol]): Poses of the molecules to be evaluated.
        scores (list[float]): Docking scores of the poses.
        target (str): Name of the target protein.
    """
    def __init__(self, poses:list[Chem.Mol], scores:list[float], target:str):
        self.poses = poses
        self.lens = len(self.poses)
        self.scores = scores
        if self.lens != len(scores):
            raise ValueError(f"Length of poses ({self.poses}) and scores ({len(scores)}) must be equal.")
        self.target = target

    def interactions(self):
        """Calculate interactions for all molecules in the target."""
        target = self.target
        key_inters = allkey_inters[target]
        empty_pdb = next((TARGET_PATH / target).glob(f'{target}*.pdb'))
        return interactions(poses=self.poses, empty_pdb=empty_pdb, key_inters=key_inters)
    
    def ligand_efficiency(self):
        """Calculate ligand efficiency for all molecules in the target."""
        le_results = process_map(calc_le, zip(self.poses, self.scores),
                                 desc=f"Calculating ligand efficiency for {self.target}...", chunksize=10)
        return le_results
    
    def sucos(self):
        """Calculate SUCOS for all molecules in the target."""
        target = self.target
        mol_true = read_sdf(next((TARGET_PATH / target).glob(f'{target}*.sdf')))
        sucos_fn = partial(check_sucos, mol_true=mol_true)
        sucos_scores = process_map(sucos_fn, self.poses,
                                   desc=f"Calculating SUCOS for {target}...", chunksize=10,)
        return sucos_scores

    def eval(self):
        """Evaluate docking results for the target."""
        # Calculate all metrics
        inters = self.interactions()
        le_results = self.ligand_efficiency()
        sucos_scores = self.sucos()
        
        eval_results = [{'pose': pose, 'score': score, **le, **sucos, **inter}
                 for pose, score, inter, le, sucos in 
                 zip(self.poses, self.scores, inters, le_results, sucos_scores)]
        return eval_results


def extract_results(results:list[dict], mode:Literal['dock', 'score_only']) -> tuple[list, list]:
    """Extract docking results from the results directory."""
    mapping = {
        'dock': ('pose', 'docking score'),
        'score_only': ('mol', 'scoring score')
    }
    p, s = mapping[mode]
    poses = [res[p] for res in results]
    scores = [res[s] for res in results]
    return poses, scores


def dock_eval(work_dir:Path, target:str) -> list[dict]:
    """Evaluate docking results for a specific target protein.
    Args:
        work_dir (Path): Path to the working directory.
        target (str): Name of the target protein.
    Returns:
        list[dict]: Evaluated results for the target.
    """

    def _read_and_validate(results_dir, target,lens, mode, error_msg):
        pkl_file = find_dockpkl(results_dir, mode=mode)
        if not pkl_file:
            raise RuntimeError(f"Results not found for {target}, {error_msg}.")
        
        results = read_pkl(pkl_file)
        if len(results) < lens:
            raise RuntimeError(f"Incomplete results for {target}, {error_msg}.")
        return results

    def _extract_and_eval(results, mode):
        poses, scores = extract_results(results, mode=mode)
        return DockEval(poses, scores, target).eval()

    # Initialize directories and read molecules
    target_dir = work_dir / target
    results_dir = target_dir / 'results'
    _, mols = read_in(target_dir)
    lens = len(mols)

    # Check before evaluation
    dock_results = _read_and_validate(
        results_dir, target, lens, mode='dock',
        error_msg="please run `tarpass dock -mode dock`."
    )
    
    if conformation_check(mols): # Check scoring results if original molecules are 3D
        score_results = _read_and_validate(
            results_dir, target, lens, mode='score_only',
            error_msg="please run `tarpass dock -mode score_only`."
        )
    else:
        scorepkl = None

    # Extract results and evaluate
    dock_eval_results = _extract_and_eval(results=dock_results, mode='dock')
    score_eval_results = _extract_and_eval(results=score_results, mode='score_only') if scorepkl else []
    
    #Combine results
    eval_all_results = []
    smis = to_smiles(mols)
    for idx, (smiles, dock_res, score_res) in enumerate(zip(smis, dock_eval_results, score_eval_results)):
        info_di = {'index': idx, 'target': target, 'smiles': smiles,}
        # Rename keys
        prefix_dock = {f"Dock {k}": v for k, v in dock_res.items()}
        prefix_score = {f"Score {k}": v for k, v in score_res.items()} if score_res else {}
        # Combine all results
        eval_all_results.append({**info_di, **prefix_dock, **prefix_score})
    project_logger.info(f"Evaluation completed for {target}.")
    return eval_all_results

############## Execution Functions ##############

def dockeval_execute(args):
    log_config(project_logger, args)
    # Preparation for reading poses
    work_dir = Path(args.path)
    for target in TARGETS:

        # Check if target directory exists
        project_logger.info(DASHLINE)
        target_dir = work_dir / target
        if not target_dir.exists():
            project_logger.warning(f"Target folder {target_dir} not found, skipping evaluation.")
            continue
        
        # Check if results exists
        project_logger.info(f'Start evaluating docking results for {target}...')
        results_dir = target_dir / 'results'
        eval_csv = results_dir / f'dock_eval_results.csv'
        if eval_csv.exists():
            project_logger.info(f"Evaluation results already exist for {target}, skipping evaluation.")
            continue
        
        # Execute evaluation and save results
        project_logger.info(DASHLINE)
        pd.DataFrame(dock_eval(work_dir, target)).to_csv(eval_csv, index=False)
        project_logger.info(f"Evaluation results saved to {eval_csv}.")
        project_logger.info(DASHLINE)