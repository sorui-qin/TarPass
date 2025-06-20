'''
Author: Rui Qin
Date: 2025-03-10 19:34:16
LastEditTime: 2025-06-20 17:21:56
Description: 
'''
import subprocess
from typing import Literal
from utils.constant import ROOT, Mol
from utils.io import read_sdf, temp_manager


class BaseDockTask():
    """
    Running docking task with AutoDock-Vina.  
    Please note that the code is designed to perform docking with **specific targets in the benchmark** for high efficiency,  
    using the small molecules and prepared files of targets as input.  
    Therefore, the code does not possess the capability for general application in docking with other targets.

    Args:
        ligand (str): Filename of prepared sdf file or a pdbqt string.
        target (str): Target name of ligand generated for.
        mode (str, optional): Docking mode ('dock' or 'score_only'). Defaults to 'dock'.
    """
    def __init__(self, ligand:str, target:str, mode:Literal['dock', 'score_only']='dock'):
        if not isinstance(ligand, str):
            raise ValueError("Invalid ligand input, please provide a filename or pdbqt string.")
        self.ligand, self.target, self.mode = ligand, target, mode
        self.target_dir = ROOT / f"Targets/{target}"

        if not self.target_dir.exists():
            raise FileNotFoundError(f"Target Fold unfound: {self.target_dir}")
        if mode not in ('dock', 'score_only'):
            raise ValueError("Invalid docking mode, choose 'dock' or 'score_only'")


class GninaDock(BaseDockTask):
    """
    Running docking task with Gnina. Ligand should be prepared in sdf format.
    Args:
        ligand (str): Filename of prepared sdf file or a pdbqt string.
        target (str): Target name of ligand generated for.
        mode (str, optional): Docking mode ('dock' or 'score_only'). Defaults to 'dock'.
    """
    def __init__(self, ligand:str, target:str, mode:Literal['dock', 'score_only']='dock'):
        super().__init__(ligand, target, mode)
        if not ligand.endswith('.sdf'):
            raise ValueError("Invalid ligand input, please provide a sdf file.")
        try:
            self.config = next((ROOT / "dock/paras").glob(f'{self.target}.txt'))
        except StopIteration:
            raise FileNotFoundError(f"No config files {self.target}.txt found in `<INSTALL_PATH>/dock/paras`")

    def _get_result(self, output_sdf:str, in_batch:bool=False):
        mols = read_sdf(output_sdf)
        if isinstance(mols, list) and in_batch:
            return [(mol, float(mol.GetProp('minimizedAffinity'))) for mol in mols]
        else:
            mol = mols[0] if isinstance(mols, list) else mols
            return mol, float(mol.GetProp('minimizedAffinity'))

    def run(self, seed=0, exhaust=8, n_poses=1, verbose=0, in_batch=False) -> tuple[Mol, float] | list[tuple[Mol, float]]:
        """Running Gnina docking task.
        Args:
            seed (int, optional): Random seed (default: 0; ramdomly choosed)
            exhaust (int, optional): Exhaustiveness of docking. Defaults to 32.
            n_poses (int, optional): Mumber of pose to generate. Defaults to 1.
            verbosity (int, optional): Verbosity. 0: not output, 1: verbose. Defaults to 1.
            in_batch (bool, optional): docking in batch mode. Defaults to False.

        Returns:
            Tuple[Chem.Mol, float], float]: Best docking pose in RdMol object and its score.
        """
        rec_pdb = next(self.target_dir.glob("*rec*.pdb"))
        # Run docking--
        with temp_manager('.sdf') as tmp_file:
            command = [
                'gnina',
                '-r', str(rec_pdb),
                '-l', str(self.ligand),
                '--config', self.config,
                '-o', tmp_file,
                '--num_modes', str(n_poses),
                '--seed', str(seed),
                '--exhaustiveness', str(exhaust),
                '--cnn', 'fast', # using fast cnn for no GPU decive
            ]
            if self.mode == 'score_only':
                command.append('--score_only')
            devnull = None if verbose != 0 else subprocess.DEVNULL
            subprocess.run(command, check=True, stdout=devnull, stderr=devnull)
            return self._get_result(tmp_file, in_batch)