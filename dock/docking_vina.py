'''
Author: Rui Qin
Date: 2025-03-01 15:57:14
LastEditTime: 2025-06-18 19:03:16
Description: 
'''
# Adapted from https://github.com/guanjq/targetdiff/blob/main/utils/evaluation/docking_vina.py
from meeko import PDBQTMolecule, RDKitMolCreate
from typing import Tuple, List
from vina import Vina
from utils.constant import Mol, ROOT
from utils.docking import sdf2centroid
from dock.docking_gnina import BaseDockTask
import json

class VinaDock(BaseDockTask):
    """
    Running docking task with AutoDock-Vina.  
    """
    def __init__(self, ligand, target, mode='dock'):
        super().__init__(ligand, target, mode)
        self.maps = ROOT / f"dock/maps/{target}/{target}"
        if not ligand.startswith('REMARK'): # Check if input ligand is a PDBQT string
            raise ValueError("Invalid ligand input, please provide an PDBQT string.")

    def _get_center(self) -> List[float]:
        """Get center of docking grid.
        """
        try: # For apo structure
            return sdf2centroid(next(self.target_dir.glob("*ligand*.sdf")))
        except: # For holo and AlphaFold structure
            return json.loads((self.target_dir/"center.json").read_text())["center"]

    def maps_prep(self, box_size:list=[20, 20, 20]):
        """Preparation for affinity map files for given target except `HDAC6`.  
           HDAC6 must be prepared with `Target/HDAC6/preprocess.ipynb`.

        Args:
            box_size (list, optional): The size of docking grid (Å). Defaults to [20, 20, 20].
        """
        if self.target == 'HDAC6':
            raise RuntimeError('Preparation for HDAC6 should follow `Target/HDAC6/preprocess.ipynb`')
        try:
            rec_pdbqt = next((ROOT / "dock/pdbqtfiles").glob(f"{self.target}_*.pdbqt"))
        except StopIteration:
            raise FileNotFoundError(f"No PDBQT files found for {self.target}") from None
        
        # Preparation of affinity maps
        (self.maps.parent).mkdir(exist_ok=True)
        v = Vina()
        v.set_receptor(str(rec_pdbqt))
        v.compute_vina_maps(self._get_center(), box_size)
        v.write_maps(map_prefix_filename=str(self.maps), overwrite=True)

    def run(self, seed=0, exhaust=8, n_poses=1, verbose=0) -> Tuple[Mol, float]:
        """Running AutoDock-Vina.

        Args:
            seed (int, optional): Random seed (default: 0; ramdomly choosed)
            exhaust (int, optional): Exhaustiveness of docking. Defaults to 8.
            n_poses (int, optional): Mumber of pose to generate. Defaults to 1.
            verbose (int, optional): Verbosity. 0: not output, 1: normal, 2: verbose (default: 0)

        Returns:
            Tuple[Chem.Mol, float], float]: Best docking pose in RdMol object and its score.
        """
        # Recprep
        if not any(self.maps.parent.glob("*.map")):
            self.maps_prep()
            if not any(self.maps.parent.glob("*.map")):
                raise RuntimeError("Failed to generate affinity maps")
        # Docking
        sf_name = 'ad4' if self.target == 'HDAC6' else 'vina'
        v = Vina(sf_name=sf_name, cpu=0, seed=seed, verbosity=verbose)
        v.set_ligand_from_string(self.ligand)
        v.load_maps(str(self.maps))
        if self.mode == 'score_only':
            pdbqt_mol = PDBQTMolecule(self.ligand, skip_typing=True)
            pose = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
            return pose, v.score()[0]
        elif self.mode == 'dock':
            v.dock(exhaustiveness=exhaust, n_poses=n_poses)
            score = v.energies(n_poses=n_poses)[0][0]
            pdbqt_mol = PDBQTMolecule(v.poses(), skip_typing=True)
            pose = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
            return pose, score
        else:
            raise ValueError("Invalid docking mode, choose 'dock' or 'score_only'")