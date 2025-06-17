'''
Author: Rui Qin
Date: 2025-03-15 13:52:13
LastEditTime: 2025-06-17 21:42:10
Description: 
'''
import argparse
from tqdm import tqdm
from utils.io import read_pkl, append_pkl, temp_manager
from utils.logger import project_logger, log_config
from utils.docking import LigPrep
from utils.constant import ROOT, TARGETS, DASHLINE
from utils.preprocess import read_in
from pathlib import Path

class Dock():
    def __init__(self, mol, target, args):
        self.mol = mol
        self.idx = mol.GetProp('_Name')
        if self.mol.GetNumConformers() == 0 and args.mode == 'score_only':
            raise ValueError(f"Conformation is not detected, `score_only` mode is banned.")
        self.target = target
        self.method_name = args.method
        if self.method_name == 'vina':
            from dock.docking_vina import VinaDock
            self.method = VinaDock
        else:
            from dock.docking_gnina import GninaDock
            self.method = GninaDock
        self.args = args
    
    def ligprep(self):
        prep = LigPrep(self.mol, self.args.reset)
        if self.method_name == 'vina':
            return prep.get_pdbqt()
        with temp_manager('.sdf', auto_remove=False) as sdf_file:
            if prep.save_sdf(sdf_file):
                return sdf_file
            else:
                return None

    def run(self):
        ligand = self.ligprep()
        if ligand is None:
            project_logger.error(f"Failed to prepare ligand {self.idx}.")
            return None, None
        try:
            dock = self.method(ligand, self.target, self.args.mode)
            return dock.run(
                seed=self.args.seed,
                exhaust=self.args.exhaust,
                n_poses=self.args.poses,
                verbose=self.args.verbose,
                )
        except Exception as e:
            project_logger.error(f"Docking failed for {self.idx}.")
            project_logger.error(e)
            return None, None
        finally:
            if self.method_name != 'vina':
                Path(ligand).unlink()

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
    return -1

############## Execution Functions ##############

def setup_arguments(parser: argparse.ArgumentParser):
    #group1 = parser.add_argument_group("Necessary arguments")
    #group1.add_argument('-p', '--path', required=True, type=str, help='path to the folder where generated molecules for testing will be stored.')

    group2 = parser.add_argument_group("Docking parameters\ndefaultly loaded from config file `configs/dock/gnina_dock.yml`")
    group2.add_argument('-m', '--method', type=str, help='docking method to use (`gnina` or `vina`).')
    group2.add_argument('--verbose', type=int, help='verbosity level of the docking process.')
    group2.add_argument('--seed', type=int, help='random seed for docking.')
    group2.add_argument('--exhaust', type=int, help='exhaustiveness of docking.')
    group2.add_argument('--poses', type=int, help='number of poses to generate.')
    group2.add_argument('--mode', type=str, choices=['dock', 'score_only'], help='docking mode.')
    group2.add_argument('--config', type=str, default=f'{ROOT}/configs/dock/gnina_dock.yml', help='path to the configuration file')

    group3 = parser.add_argument_group("Optional arguments")
    group3.add_argument('--reset', action="store_true", help="reset the original 3D conformation if available")
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
        # Breakpoint check
        Path(target_dir/'results').mkdir(parents=True, exist_ok=True)
        result_pkl = target_dir/f'results/{args.method}-{args.mode}_docking_results.pkl'
        latest_idx = breakpoint_check(result_pkl, total_lens)
        # Docking
        for i, mol in tqdm(enumerate(latest:=mols[latest_idx+1:]), 
                           desc=f'Docking with {target}', total=len(latest)):
            index = range(total_lens)[4+i]+1
            dock = Dock(mol, target, args)
            pose, score = dock.run()
            # Save results
            if args.mode == 'dock':
                append_pkl(result_pkl, 
                        {'index': index, 'mol': mol, 'pose': pose, 'docking score': score})
            elif args.mode == 'score_only':
                append_pkl(result_pkl, 
                        {'index': index, 'mol': mol, 'vina score': score})
        
        project_logger.info(f'Docking in {target} completed. Results saved in {result_pkl}.')
        print(DASHLINE)