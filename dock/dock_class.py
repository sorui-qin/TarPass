'''
Author: Rui Qin
Date: 2025-06-20 16:16:54
LastEditTime: 2025-06-23 16:32:40
Description: 
'''
from abc import ABC, abstractmethod
from utils.constant import Path, Mol
from utils.io import temp_manager, write_sdf
from utils.docking import LigPrep
from utils.logger import project_logger
from utils.preprocess import conformation_check
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')  # type: ignore

class DockBase(ABC):
    SUPPORTED_METHODS = {
        'vina': ('dock.docking_vina', 'VinaDock'),
        'gnina': ('dock.docking_gnina', 'GninaDock')
    }

    def __init__(self, target:str, args):
        self.target = target
        self.args = args
        self.method_name = args.method
        self.method = self._load_method()
    
    def _load_method(self):
        if self.method_name not in self.SUPPORTED_METHODS:
            raise ValueError(f"Unsupported method: {self.method_name}")
        
        module_path, class_name = self.SUPPORTED_METHODS[self.method_name]
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
        preps = [LigPrep(mol, self.args.reset).ligprep() for mol in self.mols]
        
        with temp_manager('.sdf', auto_remove=False) as sdf_file:
            write_sdf(sdf_file, preps) if preps else None
            if not preps:
                return None
            return sdf_file

    def _docking(self, ligand: str):
        #TODO: Add progress bar for batch docking mode
        return super()._docking(ligand, in_batch=True)