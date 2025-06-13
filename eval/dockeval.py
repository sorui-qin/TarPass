'''
Author: Rui Qin
Date: 2025-06-13 11:22:44
LastEditTime: 2025-06-13 20:55:23
Description: 
'''
import argparse
from eval.module.druglikeness import molwt
from eval.module.sucos import check_sucos
from functools import partial
from rdkit import Chem
from pathlib import Path
from tqdm.contrib.concurrent import process_map
from interaction.interaction_tools import allkey_inters, interactions
from utils.constant import ROOT, TARGETS, TMP
from utils.eval import find_dockpkl
from utils.io import read_pkl, read_sdf
from utils.logger import log_config, project_logger
from utils.preprocess import conformation_check, read_in, standard_mol
import numpy as np

TARGET_PATH = ROOT / 'Targets'

class DockEval():
    def __init__(self, poses:list[Chem.Mol], scores:list[float], target:str):
        self.poses = poses
        self.lens =len(self.poses)
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
        def calc_le(args:tuple[Chem.Mol, float]) -> dict:
            pose, score = args
            mol = standard_mol(pose)
            if score <= 0:
                le_heavyatom = mol.GetNumHeavyAtoms() / abs(score) if abs(score) > 0 else np.nan
                le_mw = molwt(mol) / abs(score) if abs(score) > 0 else np.nan
            return {'LE_heavyatom': le_heavyatom, 'LE_mw': le_mw}
        
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



def dock_eval(work_dir:Path, target:str):
    target_dir = work_dir / target
    _, mols = read_in(target_dir)

    # Check before evaluation
    if not (dockpkl:=find_dockpkl(target_dir / 'results', mode='dock')):
            raise RuntimeError(f"Docking results not found for {target}, please run `tarpass dock` first.")
    dock_results = read_pkl(dockpkl)
    if len(dock_results) < len(mols):
        raise RuntimeError(f"Incomplete docking results for {target}, please rerun `tarpass dock`.")


def execute(args):
    log_config(project_logger, args)
    # Preparation for reading poses
    work_dir = Path(args.path)
    for target in TARGETS:
        pass