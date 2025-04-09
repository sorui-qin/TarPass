'''
Author: Rui Qin
Date: 2025-03-15 13:52:13
LastEditTime: 2025-04-09 15:21:51
Description: 
'''
import argparse
from tqdm import tqdm
from utils.io import read_pkl, append_pkl, temp_manager
from utils.logger import project_logger
from utils.docking import LigPrep
from utils.constant import TARGETS, DASHLINE
from utils.preprocess import read_in
from dock.docking_vina import VinaDock
from dock.docking_gnina import GninaDock
from pathlib import Path

class Dock():
    def __init__(self, mol, target, args):
        self.mol = mol
        self.target = target
        self.method = VinaDock if args.method == 'vina' else GninaDock
        self.args = args
    
    def ligprep(self):
        prep = LigPrep(self.mol, self.args.optimize, self.args.reset)
        if self.method == VinaDock:
            return prep.get_pdbqt()
        with temp_manager('.sdf', auto_remove=False) as sdf_file:
            prep.save_sdf(sdf_file)
            return sdf_file

    def run(self):
        ligand = self.ligprep()
        try:
            dock = self.method(ligand, self.target, self.args.mode)
            return dock.run(
                seed=self.args.seed,
                exhaust=self.args.exhaust,
                n_poses=self.args.poses,
                verbose=self.args.verbose
                )
        except Exception as e:
            project_logger.error(f"Docking failed for {self.mol.GetProp('_Name')}.")
            project_logger.error(e)
            return None, None
        finally:
            if self.method == GninaDock:
                Path(ligand).unlink()

def setup_arguments(parser: argparse.ArgumentParser):
    parser.add_argument('-p', '--path', required=True, type=str, help='Path to the folder where generated molecules for testing will be stored.')
    parser.add_argument('--method', type=str, default='gnina', help='Docking method to use (`gnina` or `vina`). Default is gnina.')
    parser.add_argument('--verbose', type=int, default=0, help='Verbosity level of the docking process.')
    parser.add_argument('--seed', type=int, default=0, help='Random seed for docking.')
    parser.add_argument('--exhaust', type=int, default=8, help='Exhaustiveness of docking.')
    parser.add_argument('--poses', type=int, default=1, help='Number of poses to generate.')
    parser.add_argument('--mode', type=str, default='dock', help='Docking mode (`dock` or `score_only`). Default is dock.')
    parser.add_argument('--optimize', action="store_true", help="Optimized the original 3D conformation if available.")
    parser.add_argument('--reset', action="store_true", help="Reset the original 3D conformation if available.")
    return parser

def breakpoint_check(result_pkl: Path, total_lens: int) -> int:
    if result_pkl.exists():
        project_logger.info(f"Detected previous docking results.")
        latest_idx = read_pkl(result_pkl)[-1]['index']
        if latest_idx > total_lens:
            raise ValueError(f"Index {latest_idx} exceeds the total number of ligands {total_lens}.")
        project_logger.info(f"Docking will start from {latest_idx+1} of {total_lens} molecules.")
        return latest_idx
    return 0

def execute(args):
    work_dir = Path(args.path)
    print(DASHLINE)
    for target in TARGETS:
        target_dir = work_dir / target
        if not target_dir.exists():
            project_logger.warning(f"Target folder {target_dir} not found.")
            continue
        project_logger.info(f'Start performing docking in {target}...')

        # Preprocess
        _, mols = read_in(target_dir)
        print(DASHLINE)
        # Breakpoint check
        Path(target_dir/'results').mkdir(parents=True, exist_ok=True)
        result_pkl = target_dir/f'results/{args.method}_docking_results.pkl'
        latest_idx = breakpoint_check(result_pkl, len(mols))
        # Docking
        for mol in tqdm(mols[latest_idx:], desc=f'Docking with {target}', total=len(mols[latest_idx:])):
            dock = Dock(mol, target, args)
            pose, score = dock.run()
            # Save results
            append_pkl(result_pkl, [{'index': mols.index(mol), 'pose': pose, 'score': score}])
        project_logger.info(f'Docking in {target} completed. Results saved in {result_pkl}.')
        print(DASHLINE)