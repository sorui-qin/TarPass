from utils.dock import sdf2centroid, LigPrep
from pathlib import Path
import subprocess

class BaseDockTask():
    """
    Running docking task with AutoDock-Vina.  
    Please note that the code is designed to perform docking with **specific targets in the benchmark** for high efficiency,  
    using the small molecules and prepared affinity map files of targets as input.  
    Therefore, the code does not possess the capability for general application in docking with other targets.

    Args:
        ligand (str): Filename of a sdf file (Path) or a SMILES sequenece.
        target (str): Target name of ligand generated for.
        mode (str, optional): Docking mode ('dock' or 'score_only'). Defaults to 'dock'.
    """
    def __init__(self, ligand:str, target:str, mode='dock'):
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
        self.config = next(Path(f"dock/config").glob(f'{self.target}*.txt'))


p = Path('../Targets/BRD4-holo')
rec_pdb = '../Targets/BRD4-holo/BRD4-holo_2oss_rec_A.pdb'
lig_sdf = '../Targets/BRD4/BRD4_8pxa_ligand_A.sdf'
command = f"""gnina -r {rec_pdb} \
            -l {lig_sdf} \
            --config {p/'config.txt'} \
            -o {p/'redock-gnina.sdf'}"""