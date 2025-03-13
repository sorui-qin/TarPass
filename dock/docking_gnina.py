'''
Author: Rui Qin
Date: 2025-03-10 19:34:16
LastEditTime: 2025-03-13 21:01:52
Description: 
'''
from utils.io import read_sdf
from pathlib import Path
import subprocess
import tempfile


class BaseDockTask():
    """
    Running docking task with AutoDock-Vina.  
    Please note that the code is designed to perform docking with **specific targets in the benchmark** for high efficiency,  
    using the small molecules and prepared files of targets as input.  
    Therefore, the code does not possess the capability for general application in docking with other targets.

    Args:
        ligand (str): Filename of prepared sdf file or pdbqt string.
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
    Running docking task with Gnina.  

    """
    def __init__(self, ligand, target, mode='dock'):
        super().__init__(ligand, target, mode)
        try:
            self.config = next(Path("dock/config").glob(f'{self.target}.txt'))
        except StopIteration:
            raise FileNotFoundError(f"No config files {self.target}.txt found in `dock/config`")
    
    def _get_result(self, output_sdf:str):
        mol = read_sdf(output_sdf)
        return mol, float(mol.Getprop('minimizedAffinity'))
    
    def run(self, output_sdf:str):
        rec_pdb = next(self.target_dir.glob("*rec*.pdb"))
        with tempfile.TemporaryDirectory() as tmp:
            command = f"""gnina -r {rec_pdb} \
                        -l {self.ligand} \
                        --config {self.config} \
                        -o {tmp}/{output_sdf}"""
            subprocess.run(command, shell=True)
            return self._get_result(f"{tmp}/{output_sdf}")
        TODO