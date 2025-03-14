'''
Author: Rui Qin
Date: 2025-03-10 19:34:16
LastEditTime: 2025-03-14 20:31:23
Description: 
'''
from typing import Optional, Tuple, List
from utils.io import read_sdf, temp_manager
from pathlib import Path
import rdkit.Chem as Chem
import subprocess


class BaseDockTask():
    """
    Running docking task with AutoDock-Vina.  
    Please note that the code is designed to perform docking with **specific targets in the benchmark** for high efficiency,  
    using the small molecules and prepared files of targets as input.  
    Therefore, the code does not possess the capability for general application in docking with other targets.

    Args:
        ligand (str): preparedFilename of prepared sdf file or pdbqt string.
        target (str): Target name of ligand generated for.
        mode (str, optional): Docking mode ('dock' or 'score_only'). Defaults to 'dock'.
    """
    def __init__(self, ligand:str, target:str, mode='dock'):
        if not isinstance(ligand, str):
            raise ValueError("Invalid ligand input, please provide a filename or pdbqt string.")
        self.ligand, self.target, self.mode = ligand, target, mode
        self.target_dir = Path(f"Targets/{target}")

        if not self.target_dir.exists():
            raise FileNotFoundError(f"Target Fold unfound: {self.target_dir}")
        if mode not in ('dock', 'score_only'):
            raise ValueError("Invalid docking mode, choose 'dock' or 'score_only'")


class GninaDock(BaseDockTask):
    """
    Running docking task with Gnina. Ligand should be prepared in sdf format.
    """
    def __init__(self, ligand, target, mode='dock'):
        super().__init__(ligand, target, mode)
        if not ligand.endswith('.sdf'):
            raise ValueError("Invalid ligand input, please provide a sdf file.")
        try:
            self.config = next(Path("dock/config").glob(f'{self.target}.txt'))
        except StopIteration:
            raise FileNotFoundError(f"No config files {self.target}.txt found in `dock/config`")


    def _get_result(self, output_sdf:str):
        mol = read_sdf(output_sdf)
        return mol, float(mol.GetProp('minimizedAffinity'))

    def run(self, seed=0, exhaust=8, n_poses=1, verbose=0) -> Tuple[Optional[List[Chem.Mol]], float]:
        """Running Gnina docking task.
        Args:
            seed (int, optional): Random seed (default: 0; ramdomly choosed)
            exhaust (int, optional): Exhaustiveness of docking. Defaults to 32.
            n_poses (int, optional): Mumber of pose to generate. Defaults to 1.
            Verbosity (int, optional): Verbosity. 0: not output, 1: verbose. Defaults to 1.

        Returns:
            Tuple[Optional[List[Chem.Mol]], float]: Docking score and poses list in RdMol object.
        """
        output_dir = './tmp'
        Path(output_dir).mkdir(exist_ok=True)
        rec_pdb = next(self.target_dir.glob("*rec*.pdb"))
        # Run docking
        with temp_manager('.sdf', output_dir) as tmp_file:
            command = [
                'gnina',
                '-r', rec_pdb,
                '-l', self.ligand,
                '--config', self.config,
                '-o', tmp_file,
                '--num_modes', str(n_poses),
                '--seed', str(seed),
                '--exhaustiveness', str(exhaust),
            ]
            if self.mode == 'score_only':
                command.append('--score')
            if not verbose:
                command.append('--quiet')
            subprocess.run(command, check=True)
            return self._get_result(tmp_file)

if __name__ == '__main__':
    ligand = 'testfile/5ht2a_test.sdf'
    target = '5HT2A'
    dock = GninaDock(ligand, target)
    mol, score = dock.run()
    print(score)