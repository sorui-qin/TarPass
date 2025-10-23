'''
Author: Rui Qin
Date: 2025-06-20 16:16:54
LastEditTime: 2025-10-23 14:27:45
Description: 
'''
import argparse
from abc import ABC, abstractmethod
from utils.constant import Path, Mol
from utils.io import temp_manager, write_sdf
from utils.docking import LigPrep
from utils.logger import project_logger
from utils.preprocess import conformation_check
from rdkit import RDLogger
from tqdm.contrib.concurrent import process_map
RDLogger.DisableLog('rdApp.*')  # type: ignore

METHODS = {
    'vina': ('dock.docking_vina', 'VinaDock'),
    'gnina': ('dock.docking_gnina', 'GninaDock')
}

class DockBase(ABC):
    def __init__(self, target:str, args):
        self.target = target
        self.args = args
        self.method_name = args.method
        self.method = self._load_method()
    
    def _load_method(self):
        if self.method_name not in METHODS:
            raise ValueError(f"Unsupported method: {self.method_name}")
        
        module_path, class_name = METHODS[self.method_name]
        module = __import__(module_path, fromlist=[class_name])
        return getattr(module, class_name)

    @abstractmethod
    def _ligprep(self):
        pass

    def _docking(self, ligand: str, **kwargs):
        dock = self.method(ligand, self.target, self.args.mode)
        try:
            return dock.run(
                seed=self.args.seed,
                exhaust=self.args.exhaust,
                n_poses=self.args.poses,
                verbose=self.args.verbose,
                **kwargs
            )
        finally:
            if self.method_name != 'vina':
                Path(ligand).unlink()

    def run(self):
        """Running the docking task.
        """
        ligand = self._ligprep()
        if ligand is None:
            raise RuntimeError("Failed to prepare ligand.")
        else:
            return self._docking(ligand)


class SingleDock(DockBase):
    def __init__(self, mol:Mol, target:str, args):
        super().__init__(target, args)
        self.mol = mol
        self.idx = mol.GetProp('_Name')
        if self.mol.GetNumConformers() == 0 and args.mode == 'score_only':
            raise ValueError(f"Conformation is not detected, `score_only` mode is banned.")
    
    def _ligprep(self) -> str | None:
        prep = LigPrep(self.mol, self.args.reset)
        
        if self.method_name == 'vina':
            return prep.get_pdbqt()
        with temp_manager('.sdf', auto_remove=False) as sdf_file:
            return sdf_file if prep.save_sdf(sdf_file) else None

def _prep(args):
    (index, mol), reset = args
    # NOTE: Reset the name of the ligand because it will be lost in pickle procession during multiprocessing.
    mol.SetProp('_Name', str(index))
    return LigPrep(mol, reset).ligprep()

class BatchDock(DockBase):
    def __init__(self, mols:list[Mol], target:str, args):
        super().__init__(target, args)
        self.mols = mols
        if not conformation_check(self.mols) and args.mode == 'score_only':
            raise ValueError(f"Conformation is not detected, `score_only` mode is banned.")
        if self.method_name == 'vina':
            raise ValueError(f"Batch docking with Vina is not supported.")

    def _ligprep(self) -> str | None:
        project_logger.info(f"Preparing ligands for {self.target}...")
        mol_tuples = [((index, mol), self.args.reset) for index, mol in enumerate(self.mols)]
        preps = process_map(_prep, mol_tuples, desc='Executing LigPrep')
        
        with temp_manager('.sdf', auto_remove=False) as sdf_file:
            write_sdf(sdf_file, preps) if preps else None
            if not preps:
                return None
            return sdf_file

    def _docking(self, ligand: str):
        #TODO: Add progress bar for batch docking mode
        return super()._docking(ligand, in_batch=True)


def easydock(mols:list[Mol], target:str):
    """A simple interface for docking using default parameters with Gnina.
    Args:
        mols (list[Chem.Mol]): List of ligands to be docked.
        target (str): Target name.
    """
    args = argparse.Namespace(
            method='gnina',
            mode='dock',
            seed=0,
            exhaust=8,
            poses=1,
            verbose=False,
            reset=False
        )
    dock = BatchDock(mols, target, args)
    return dock.run()