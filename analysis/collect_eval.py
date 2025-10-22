'''
Author: Rui Qin
Date: 2025-07-07 17:25:29
LastEditTime: 2025-10-21 21:15:34
Description: 
'''
from dataclasses import dataclass, field
from typing import Any, Literal, Optional

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from utils.constant import TARGETS, Mol, Path
from utils.eval import find_dockpkl, read_eval
from utils.io import read_pkl, write_pkl
from utils.logger import project_logger
from utils.preprocess import conformation_check, to_mols


class CollectEval:
    """
    Collect raw evaluation results for docking and molecule evaluation.
    """
    def __init__(self, results_dir: str|Path):
        self.dir = Path(results_dir)
        self.target = self.dir.parent.name
        self.mols = self._get_mols()
        self.confs = conformation_check(self.mols)
        self.dock_eval, self.mole_eval = self._read_eval()
        self.dock_dicts = [res['dock'] for res in self.dock_eval]
        self.score_dicts = [res['score'] for res in self.dock_eval] if self.confs else []
    
    def _get_mols(self) -> list:
        dockpkl = find_dockpkl(self.dir, mode='dock')
        return [res['mol'] for res in read_pkl(dockpkl)]
    
    def _read_eval(self) -> tuple[list[dict], pd.DataFrame]:
        dock_eval = read_eval(self.dir, mode='dock')
        mole_eval = read_eval(self.dir, mode='mole')
        if len(dock_eval) != len(mole_eval):
            raise ValueError(f"Docking and Molecule evaluation results length mismatch for {self.target}.")
        return dock_eval, mole_eval
    
    def split_dock(self, mode:Literal['dock', 'score'], 
                    split_key:str='matched_rate') -> tuple[dict, dict]:
        """Split docking evaluation results into numerical and interaction parts.

        Args:
            mode (str, optional): Reading mode. Defaults to Literal['dock', 'score'].
            split_key (str, optional): Key to split. Defaults to 'matched_rate'.
        """
        dicts = self.dock_dicts if mode == 'dock' else self.score_dicts
        if not dicts:
            return {}, {}
        
        if split_key not in dicts[0]:
            raise ValueError(f"Split key '{split_key}' not found in dictionary. Available keys: {list(dicts[0].keys())}")
    
        keys = list(dicts[0].keys())
        try:
            split_idx = keys.index(split_key) + 1
        except ValueError:
            raise ValueError(f"Split key '{split_key}' not found in dictionary.")
        
        numericals = [{k: d[k] for k in keys[:split_idx]} for d in dicts]
        interactions = [{k: d[k] for k in keys[split_idx:]} for d in dicts]
        
        return (pd.DataFrame(numericals).to_dict(orient='list'),
                pd.DataFrame(interactions).to_dict(orient='list'))
    
    def split_mole(self) -> tuple[dict, dict, dict]:
        """Split molecule evaluation results into descriptor, structure, and alert parts.
        """
        DESCI_IDX = [*range(2, 12), *range(24, 29)] 
        STRUC_IDX = [*range(12, 18), *range(24, 36)]
        ALERT_IDX = list(range(18, 24))

        if (cols := self.mole_eval.shape[1]) != 36:
            raise ValueError(f"Unexpected number of columns: {cols}. Expected 36.")
    
        return (
        self.mole_eval.iloc[:, DESCI_IDX].to_dict(orient='list'),
        self.mole_eval.iloc[:, STRUC_IDX].to_dict(orient='list'),
        self.mole_eval.iloc[:, ALERT_IDX].to_dict(orient='list')
        )


@dataclass(frozen=True)
class MolInfo:
    """Molecular information including SMILES, InChI, RDKit Mol object, and pose."""
    rdmol: list[Mol]
    pose: list[Mol]
    smiles: list[str]
    inchikey: list[str]

@dataclass(frozen=True)
class StructInfo:
    """Structure-based information for docking or scoring.
    """
    numericals: dict[str, list] = field(default_factory=dict)
    interactions: dict[str, list] = field(default_factory=dict)

@dataclass(frozen=True)
class PropInfo:
    """Molecular properties information including descriptors, structural features, and alerts.
    """
    descriptors: dict[str, list[float]] = field(default_factory=dict)
    structural: dict[str, list[float]] = field(default_factory=dict)
    alerts: dict[str, list[bool]] = field(default_factory=dict)

@dataclass(frozen=True)
class MoleculesData:
    """
    Molecules data class to hold docking (scoring) and property information for a list of molecules.
    """
    Mol: MolInfo
    Dock: StructInfo
    Score: Optional[StructInfo]
    Prop: PropInfo

    def __getitem__(self, index: int) -> dict[str, Any]:
        """Get all detailed molecule information by index.
        """
        num_molecules = self.__len__()
        if index < 0 or index >= num_molecules:
            raise IndexError(f"Index out of range: {index}. Valid range is 0 to {num_molecules - 1}.")
        return {
            "Mol": {
                "rdmol": self.Mol.rdmol[index],
                "pose": self.Mol.pose[index],
                "smiles": self.Mol.smiles[index],
                "inchikey": self.Mol.inchikey[index]
            },
            "Dock": {
                "numerical": {k: v[index] for k, v in self.Dock.numericals.items()},
                "interactions": {k: v[index] for k, v in self.Dock.interactions.items()}
            },
            "Score": {
                "numerical": {k: v[index] for k, v in self.Score.numericals.items()},
                "interactions": {k: v[index] for k, v in self.Score.interactions.items()}
            } if self.Score else None,
            "Prop": {
                "Descriptors": {k: v[index] for k, v in self.Prop.descriptors.items()},
                "Structural": {k: v[index] for k, v in self.Prop.structural.items()},
                "Alerts": {k: v[index] for k, v in self.Prop.alerts.items()}
            }
        }
    
    def __len__(self) -> int:
        return len(self.Mol.rdmol)
    
    def __post_init__(self) -> None:
        """Validate the lengths of all fields in the data class.
        """
        num_molecules = len(self.Mol.rdmol)
        for field in self.Dock.__dict__.values():
            if isinstance(field, dict):
                for values in field.values():
                    if len(values) != num_molecules:
                        raise ValueError("Unexpected data length in Dock")
        
        if self.Score:
            for field in self.Score.__dict__.values():
                if isinstance(field, dict):
                    for values in field.values():
                        if len(values) != num_molecules:
                            raise ValueError("Unexpected data length in Score")
        
        for category in self.Prop.__dict__.values():
            for values in category.values():
                if len(values) != num_molecules:
                    raise ValueError("Unexpected data length in Prop")


def collect_eval(results_dir: Path) -> MoleculesData:
    """Collect evaluation results from the specified directory for further analysis.
    
    Args:
        results_dir (Path): Path to the results directory.
    
    Returns:
        MoleculesData: A data class containing all collected information.
    """
    if not results_dir.exists():
        FileNotFoundError(f"Results directory {results_dir} does not exist.")
    
    mol_records = read_pkl(find_dockpkl(results_dir, mode='dock'))
    rdmol_list = []
    pose_list = []
    smiles_list = []
    inchikey_list = []
    for res in mol_records:
        mol = res['mol']
        rdmol_list.append(mol)
        pose_list.append(res['pose'])
        smiles_list.append(Chem.MolToSmiles(mol))
        inchikey_list.append(Chem.MolToInchiKey(mol))
    Mol = MolInfo(
        rdmol=rdmol_list,
        pose=pose_list,
        smiles=smiles_list,
        inchikey=inchikey_list
    )

    collect = CollectEval(results_dir)

    nums, inters = collect.split_dock('dock')
    Dock = StructInfo(
        numericals=nums,
        interactions=inters
    )

    nums_s, inters_s = collect.split_dock('score', 'matched_rate')
    Score = StructInfo(
        numericals=nums_s, 
        interactions=inters_s) \
    if all([nums_s, inters_s]) else None

    descriptors, structural, alerts = collect.split_mole()
    Prop = PropInfo(
        descriptors=descriptors,
        structural=structural,
        alerts=alerts
    )

    return MoleculesData(
        Mol=Mol,
        Dock=Dock,
        Score=Score,
        Prop=Prop
        )

def collect_readable(reselts_dir: Path) -> list[dict[str, Any]]:
    """Collect evaluation results from the specified directory and convert them into a readable list of dictionaries, 
    each containing molecular information presented in molecular order.

    Args:
        reselts_dir (Path): Path to the results directory.
    """
    results = collect_eval(reselts_dir)
    return [results[i] for i in range(len(results))]


def collect_eval_all(work_dir:str|Path, prefix:str, save_dir:Optional[str|Path]=None, readable:bool=False) -> None:
    """Collect evaluation results for all targets to a pickle file.
    Args:
        work_dir (str|Path): Path to the working directory containing target folders.
        prefix (str): Prefix for the output pickle file.
        save_dir (Optional[str|Path]): Directory to save the results. If None, saves in work_dir.
    """
    work_dir = Path(work_dir)
    if not work_dir.exists():
        raise FileNotFoundError(f"Working directory {work_dir} does not exist.")
    
    results = {}
    for target in tqdm(TARGETS, desc='Collecting eval results from targets'):
        target_dir = work_dir / target
        if not target_dir.exists():
            project_logger.warning(f"Target folder {target_dir} not found, skipping evaluation.")
            continue
        results[target] = collect_eval(target_dir / 'results') if not readable else collect_readable(target_dir / 'results')
    
    if save_dir:
        Path(save_dir).mkdir(parents=True, exist_ok=True)
    else:
        save_dir = work_dir
        save_dir.mkdir(exist_ok=True)
    pkl_path = Path(save_dir) / f'{prefix}_eval_results.pkl'

    write_pkl(pkl_path, results)
    project_logger.info(f"Evaluation results collected and saved to {pkl_path}.")


class DataCroupier: # Just a guilty pleasure to use this name.
    """Collections for evaluation results, including test data, reference, and random.
    """
    def __init__(
            self,
            test_data: MoleculesData,
            reference: MoleculesData,
            rand: MoleculesData,
            target: Optional[str]=None
            ):
        self.test_data = test_data
        self.reference = reference
        self.rand = rand
        self.target = target


class AnalysisBase:
    """Base class for analysis, providing a common interface for different analyses.
    """
    def __init__(self, analysis:DataCroupier, attr_name):
        self.test = getattr(analysis.test_data, attr_name)
        self.ref = getattr(analysis.reference, attr_name)
        self.rand = getattr(analysis.rand, attr_name)
        self.croupier_keys = ['test', 'ref', 'rand']
        self.target = analysis.target

    def _split_attr(self, attr):
        """Split the specified attribute.
        """
        if not hasattr(self, attr):
            raise AttributeError(f"Attribute '{attr}' not found in {self.__class__.__name__}.")
        
        data = getattr(self, attr)
        if isinstance(data, StructInfo):
            return data.numericals, data.interactions
        elif isinstance(data, PropInfo):
            return data.descriptors, data.structural, data.alerts
        elif isinstance(data, MolInfo):
            return data.smiles # Only return SMILES list for further processing.
        else:
            raise TypeError(f"Unsupported attribute type: {type(data)}.")