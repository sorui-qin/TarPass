# Adapted from https://github.com/guanjq/targetdiff/blob/main/utils/evaluation/docking_vina.py

from meeko import PDBQTMolecule, RDKitMolCreate
from typing import Optional, Tuple, List
from vina import Vina
import json
import rdkit.Chem as Chem
from pathlib import Path
from utils.dock import sdf2centroid, LigPrep
from dock.docking_gnina import BaseDockTask


class VinaDock(BaseDockTask):
    """
    Running docking task with AutoDock-Vina.  
    Please note that the code is designed to perform docking with **specific targets in the benchmark** for high efficiency,  
    using the small molecules and prepared affinity map files of targets as input.  
    Therefore, the code does not possess the capability for general application in docking with other targets.
    """
    def __init__(self, ligand, target, mode='dock'):
        super().__init__(ligand, target, mode)
        self.maps = Path(f"dock/maps/{target}/{target}")

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
            rec_pdbqt = next(Path("dock/pdbqtfiles").glob(f"{self.target}_*.pdbqt"))
        except StopIteration:
            raise FileNotFoundError(f"No PDBQT files found for {self.target}") from None
        
        # Preparation of affinity maps
        (self.maps.parent).mkdir(exist_ok=True)
        v = Vina()
        v.set_receptor(str(rec_pdbqt))
        v.compute_vina_maps(self._get_center(), box_size)
        v.write_maps(map_prefix_filename=self.maps, overwrite=True)


    def dock(self, pdbqt_string, sf_name:str, seed=0, exhaust=32, n_poses=1, verbose=0):
        """Running docking process with ligand PDBQT string and receptor affinity map files.

        Args:
            pdbqt_string (str): PDBQT string of ligand.
            sf_name (str): Scoring function name to use (vina, or ad4). `ad4` is only used for zinc metalloprotein.
            seed (int, optional): Random seed (default: 0; ramdomly choosed)
            exhaust (int, optional): Exhaustiveness of docking. Defaults to 32.
            n_poses (int, optional): Mumber of pose to generate. Defaults to 1.
            verbose (int, optional): verbosity 0: not output, 1: normal, 2: verbose (default: 0)

        Returns:
            score, poses: Docking score and poses list in RdMol object.
        """
        v = Vina(sf_name=sf_name, cpu=0, seed=seed, verbosity=verbose)
        v.set_ligand_from_string(pdbqt_string)
        v.load_maps(str(self.maps))
        if self.mode == 'score_only':
            return v.score()[0], []
        elif self.mode == 'dock':
            v.dock(exhaustiveness=exhaust, n_poses=n_poses)
            score = v.energies(n_poses=n_poses)[0][0]
            pdbqt_mol = PDBQTMolecule(v.poses(), skip_typing=True)
            poses = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
            return score, poses
        

    def run(self, optimize=False, **kwargs) -> Tuple[float, Optional[List[Chem.Mol]]]:
        """Running AutoDock-Vina

        Args:
            optimize (bool, optional): Optimize the conformation (default: False)

        Returns:
            Tuple[float, Optional[List[Chem.Mol]]]: Docking score and poses list in RdMol object.
        """
        # Ligprep
        pdbqt_string = LigPrep(self.ligand, optimize).get_pdbqt()
        # Recprep
        if not any(self.maps.parent.glob("*.map")):
            self.maps_prep()
            if not any(self.maps.parent.glob("*.map")):
                raise RuntimeError("Failed to generate affinity maps")
        # Docking
        sf_name = 'ad4' if self.target == 'HDAC6' else 'vina'
        return self.dock(pdbqt_string, sf_name=sf_name, **kwargs)

if __name__ == '__main__':
    v = VinaDock('Targets/BRD4/BRD4_8pxa_ligand_A.sdf', 'BRD4-holo')
    score, pose = v.run(verbose=1)
    with Chem.SDWriter('Targets/BRD4-holo/redock.sdf') as w:
        for mol in pose:
            w.write(mol)
    #p = LigPrep('./Targets/BRD4/BRD4_8pxa_ligand_A.sdf')
    #print(p.get_pdbqt())