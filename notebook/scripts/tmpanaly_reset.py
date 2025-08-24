import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from module.chemdist import get_median_iqr
from module.significance import SignificanceTester
from utils.io import load_json, read_pkl

from analysis.analysis import ROOT, TARGETS, process_pli, tqdm
from analysis.collect_eval import MoleculesData, StructInfo

DATA_DIR = ROOT / 'data'
REF_PKL = DATA_DIR / 'Reference_eval_results.pkl'
DEC_PKL = DATA_DIR / 'Decoy_eval_results.pkl'

class CollectEval:
    def __init__(self, results_dir: str|Path):
        self.dir = Path(results_dir)
        self.target = self.dir.parent.name
        self.dock_eval = load_json(self.dir / 'reset_eval_results.json')
        self.dock_dicts = [res['reset'] for res in self.dock_eval]
    
    def split_reset(self,split_key:str='matched_rate') -> tuple[dict, dict]:
        dicts = self.dock_dicts
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


@dataclass(frozen=True)
class TmpData:
    Dock: StructInfo

    def get_molecule(self, index: int):
        return {
            "Dock": {
                "numerical": {k: v[index] for k, v in self.Dock.numericals.items()},
                "interactions": {k: v[index] for k, v in self.Dock.interactions.items()}
                },
            }


def collect_eval(results_dir: Path) -> TmpData:
    if not results_dir.exists():
        FileNotFoundError(f"Results directory {results_dir} does not exist.")
    collect = CollectEval(results_dir)
    nums, inters = collect.split_reset()
    Dock = StructInfo(numericals=nums, interactions=inters)
    return TmpData(Dock=Dock)


class DataCroupier:
    def __init__(
            self,
            test_data: TmpData,
            reference: MoleculesData,
            decoy: MoleculesData,
            target: Optional[str]=None
            ):
        self.test_data = test_data
        self.reference = reference
        self.decoy = decoy
        self.target = target

class AnalysisBase:
    def __init__(self, analysis:DataCroupier, attr_name):
        self.test = getattr(analysis.test_data, attr_name)
        self.ref = getattr(analysis.reference, attr_name)
        self.decoy = getattr(analysis.decoy, attr_name)
        self.croupier_keys = ['test', 'ref', 'decoy']
        self.target = analysis.target

    def _split_attr(self, attr):
        if not hasattr(self, attr):
            raise AttributeError(f"Attribute '{attr}' not found in {self.__class__.__name__}.")
        
        data = getattr(self, attr)
        if isinstance(data, StructInfo):
            return data.numericals, data.interactions
        

class PLIAnalysis(AnalysisBase):
    """Analysis for protein-ligand interactions (PLI) related metrics from docking result."""
    def __init__(self, collections: DataCroupier):
        super().__init__(collections, 'Dock')
        
        self.numericals, self.interactions = [], []
        for attr in self.croupier_keys:
            num, inter = self._split_attr(attr) # type: ignore
            self.numericals.append(num)
            self.interactions.append(inter)
        self.test_numerical = self.numericals[0]
        self.test_interaction = self.interactions[0]
    
    def _get_median_iqr(self, key: str) -> pd.DataFrame:
        med, q1, q3 = get_median_iqr(self.test_numerical[key])
        return pd.DataFrame(
            {key: f'[{round(med, 4)}, {round(q1, 4)}, {round(q3, 4)}]'},
            index=[0])
    
    def _get_mean(self, key:str) -> pd.DataFrame:
        mean = np.nanmean(self.test_numerical[key])
        return pd.DataFrame({key: round(mean, 4)}, index=[0])
    
    def _score_significance(self, key='score') -> pd.DataFrame:
        control_data = self.test_numerical[key]
        data_groups = [i[key] for i in self.numericals[1:]]
        tester = SignificanceTester(data_groups, control_data, key)
        return tester.ref_decoy_analysis(alternative='less')
    
    def _count_interactions(self):
        total = len(self.test_interaction['detected_interactions'])
        all_inters = []
        for idx in range(total):
            interaction = self.test_interaction['detected_interactions'][idx]
            if interaction is None:
                all_inters.append(0)
            for v in interaction.values():
                all_inters.append(sum(len(i) for i in v))
        
        return pd.DataFrame({
            'average_interactions': round(np.mean(all_inters), 4),
        }, index=[0])


    def analysis(self) -> pd.DataFrame:
        affinity_cols = ['score']
        boolean_cols = [
            'sucos', 'centroid shift','shape_similarity',
            'electrostatic_similarity', 'LE_heavyatom', 'LE_mw',
            'no_clashes', 'fully_matched', 'matched_rate'
            ]
        
        dfs = (
        [self._get_median_iqr(col) for col in affinity_cols] +
        [self._get_mean(col) for col in boolean_cols] +
        [self._count_interactions()] +
        [self._score_significance()]
        )
        df = pd.concat(dfs, axis=1)
        df.index = [self.target]  # type: ignore
        return df

def _load_test(path: Path):
    if path.is_file() and path.suffix == '.pkl':
        return read_pkl(path), 'pkl'
    if path.is_dir():
        return {}, 'dir'
    raise ValueError(f"Invalid path: {path}. It should be directory or a pickle file.")
    
def execute(path):
    work_dir = Path(path)
    tests_data, read_source = _load_test(work_dir)
    
    refs = read_pkl(REF_PKL)
    decoys = read_pkl(DEC_PKL)

    results = {'pli_dock': []}
    
    for target in tqdm(TARGETS, desc='Analysis targets'):
        # Extract the test data for the target
        if read_source == 'pkl':
            test_target = tests_data.get(target, None) # type: ignore
            if test_target is None:
                continue
        elif read_source == 'dir':
            if not (target_dir := work_dir / target).exists():
                continue
            test_target = collect_eval(target_dir / 'results')

        croupier = DataCroupier(test_target, refs[target], decoys[target], target) # type: ignore
        pli_info = PLIAnalysis(croupier).analysis()
        results['pli_dock'].append(pli_info)
    results_concat = {k: pd.concat(v, ignore_index=False) for k, v in results.items()}

    # PLI Analysis
    results_concat['pli_dock'] = process_pli(results_concat['pli_dock'])

    output_dir = work_dir

    for _, v in results_concat.items():
        if not v.empty:
            v.to_csv(output_dir / f'tarpass_anlysis_reset.csv', index=True)

path = sys.argv[1]
execute(path)