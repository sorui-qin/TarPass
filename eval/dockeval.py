'''
Author: Rui Qin
Date: 2025-06-13 11:22:44
LastEditTime: 2025-06-26 14:14:15
Description: 
'''
import itertools
from functools import partial
from pathlib import Path
from typing import Literal
from rdkit.Chem import Mol
from tqdm.contrib.concurrent import process_map
from interaction.interaction_tools import allkey_inters, interactions
from module.druglikeness import calc_le
from module.intermolecular_distance import (centriod_distance,
                                            check_intermolecular_distance)
from module.sucos import check_sucos
from utils.constant import DASHLINE, TARGET_PATH, TARGETS
from utils.eval import find_dockpkl
from utils.io import dump_json, read_pdb_rdmol, read_pkl, read_sdf
from utils.logger import log_config, project_logger
from utils.preprocess import conformation_check, read_in, to_smiles


class DockEval:
    """Evaluate docking results of generated molecules for a specific target.
    Args:
        poses (list[Chem.Mol]): Poses of the molecules to be evaluated.
        scores (list[float]): Docking scores of the poses.
        target (str): Name of the target protein.
    """
    def __init__(self, poses:list[Mol], scores:list[float], target:str):
        self.lens = len(poses)
        if self.lens != len(scores):
            raise ValueError(f"Length of poses ({len(poses)}) and scores ({len(scores)}) must be equal.")
        self.poses = poses
        self.scores = scores
        self.target = target
        self.pdb_path, self.lig_path = self._get_target()
        self.pdb_rdmol, self.lig_mol = read_pdb_rdmol(self.pdb_path), read_sdf(self.lig_path)
    
    def _get_target(self):
        target_ditr = TARGET_PATH / self.target
        pdb_files = list((target_ditr).glob(f'{self.target}*.pdb'))
        if not pdb_files:
            raise FileNotFoundError(f"No pdb file found for target {self.target}")
        sdf_files = list((TARGET_PATH / self.target).glob(f'{self.target}*.sdf'))
        if not sdf_files:
            raise FileNotFoundError(f"No sdf file found for target {self.target}")
        return pdb_files[0], sdf_files[0]

    def clash(self):
        """Check for clashes between the poses and the target protein."""
        clash_fn = partial(check_intermolecular_distance, mol_cond=self.pdb_rdmol)
        return process_map(
            clash_fn, self.poses, chunksize=100,
            desc=f"Checking clashes",
            )
    
    def centriod_shift(self):
        """Calculate the centroid shift of the poses relative to the ref ligand."""
        shift_fn = partial(centriod_distance, mol_ref=self.lig_mol)
        return process_map(
            shift_fn, self.poses, chunksize=100,
            desc=f"Calculating centroid shift",
            )

    def interactions(self):
        """Calculate interactions for all molecules in the target."""
        return interactions(
            poses=self.poses, empty_pdb=self.pdb_path, 
            key_inters=allkey_inters[self.target]
            )

    def ligand_efficiency(self):
        """Calculate ligand efficiency for all molecules in the target."""
        return process_map(
            calc_le, zip(self.poses, self.scores), chunksize=100,
            desc=f"Calculating ligand efficiency",
            total=self.lens
            )

    def sucos(self):
        """Calculate SUCOS for all molecules in the target."""
        sucos_fn = partial(check_sucos, mol_true=self.lig_mol)
        return process_map(
            sucos_fn, self.poses, chunksize=100,
            desc=f"Calculating SUCOS",
            )
 
    def evaluate(self):
        """Evaluate docking results for the target."""
        # Calculate all metrics
        clashes = self.clash()
        inters_results = self.interactions()
        le_results = self.ligand_efficiency()
        sucos_scores = self.sucos()
        
        # Combine all results
        eval_results = []
        for i, score in enumerate(self.scores):
            result = {
                'score': score,
                **clashes[i],
                **sucos_scores[i],
                **le_results[i],
                **inters_results[i]
            }
            eval_results.append(result)
        return eval_results


def extract_results(results:list[dict], mode:Literal['dock', 'score_only']) -> tuple[list, list]:
    """Extract docking results from the results directory."""
    mapping = {
        'dock': ('pose', 'docking score'),
        'score_only': ('mol', 'vina score')
    }
    p, s = mapping[mode]
    poses = [res[p] for res in results]
    scores = [res[s] for res in results]
    return poses, scores


def dock_eval(mols:list[Mol], target_dir:Path) -> list[dict]:
    """Evaluate docking results for a specific target protein.
    Args:
        mols (list[Mol]): List of molecules to be evaluated.
        target_dir (Path): Path to the target directory.

    Returns:
        list[dict]: Evaluated results for the target.
    """

    def _read_and_validate(results_dir, target, lens, mode, error_msg):
        pkl_file = find_dockpkl(results_dir, mode=mode)
        if not pkl_file:
            raise RuntimeError(f"Results not found for {target}, {error_msg}.")
        
        results = read_pkl(pkl_file)
        if len(results) < lens:
            raise RuntimeError(f"Incomplete results for {target}, {error_msg}.")
        return results

    def _extract_and_eval(results, mode):
        poses, scores = extract_results(results, mode=mode)
        return DockEval(poses, scores, target).evaluate()

    # Initialize directories and read molecules
    results_dir = target_dir / 'results'
    lens = len(mols)
    target = target_dir.name

    # Check before evaluation
    dock_results = _read_and_validate(
        results_dir, target, lens, mode='dock',
        error_msg="please run `tarpass dock -mode dock`."
    )
    project_logger.info(f'Now evaluating docking results for {target}...')

    if conformation_check(mols): # Check scoring results if original molecules are 3D
        score_results = _read_and_validate(
            results_dir, target, lens, mode='score_only',
            error_msg="please run `tarpass dock -mode score_only`."
        )
    else:
        score_results = None

    # Extract results and evaluate
    dock_eval_results = _extract_and_eval(results=dock_results, mode='dock')
    if score_results:
        project_logger.info('Initial molecules are 3D, checking scoring results of initial conformation...')
        score_eval_results = _extract_and_eval(results=score_results, mode='score_only')
    else:
        project_logger.info(f'Initial molecules are not 3D, skipping scoring results for {target}.')
        score_eval_results = []
    
    #Combine results
    eval_all_results = []
    smis = to_smiles(mols)
    for idx, (smiles, dock_res, score_res) in enumerate(
        itertools.zip_longest(smis, dock_eval_results, score_eval_results, fillvalue=None)
    ):
        info_di = {'index': idx, 'target': target, 'smiles': smiles,}
        info_di.update({'dock': dock_res})
        if score_res:
            info_di.update({'score': score_res})
        eval_all_results.append(info_di)
    project_logger.info(f"Evaluation completed for {target}.")
    return eval_all_results

############## Execution Functions ##############

def dockeval_execute(args):
    log_config(project_logger, args)

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
        eval_output = results_dir / f'dock_eval_results.json'
        if eval_output.exists():
            project_logger.info(f"Evaluation results already exist for {target}, skipping evaluation.")
            continue
        
        # Execute evaluation and save results
        project_logger.info(DASHLINE)
        _, mols = read_in(target_dir)
        result = dock_eval(mols, target_dir)
        dump_json(eval_output, result)
        project_logger.info(f"Evaluation results saved to {eval_output}.")
        project_logger.info(DASHLINE)