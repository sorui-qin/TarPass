'''
Author: Rui Qin
Date: 2025-07-20 16:31:27
LastEditTime: 2025-07-22 19:18:26
Description: 
'''
import numpy as np
import pandas as pd
from analysis.collect_eval import AnalysisBase, DataCroupier
from module.chemdist import get_median_iqr, ks_distance, wasserstein_ref_decoy
from module.significance import SignificanceTester
from module.similarity import Similarity

##### Analysis for molecular distance #####

class MolDisAnalysis(AnalysisBase):
    """Analysis for molecular distance."""
    def __init__(self, collections: DataCroupier):
        super().__init__(collections, 'Mol')
        self.smis_croupier = [self._split_attr(attr) for attr in self.croupier_keys]
        self.test_smis = self.smis_croupier[0]
    
    def _internal_distance(self) -> pd.DataFrame:
        simi = Similarity(self.test_smis) # type: ignore
        return pd.DataFrame([{
            'IntDiv': round(simi.intdiv(), 4),
            'IntDiv-Scaffold': round(simi.intdiv(cal_scaffold=True), 4),
            '#Circle': round(simi.circle(), 4),
        }])
        
    def _compare_distance(self, idx:int) -> pd.DataFrame:
        key_map = {
            0: 'Ref',
            1: 'Decoy'
        }
        key = key_map[idx]
        simi = Similarity(self.test_smis, self.smis_croupier[idx+1]) # type: ignore
        return pd.DataFrame([{
            f'{key}_similarity': round(simi.similarity(), 4),
            f'{key}_scaff_similarity': round(simi.scaffold_similarity(), 4),
            f'{key}_FCD': round(simi.fcd(), 4)
        }])
    
    def analysis(self):
        dfs = (
            [self._internal_distance()] +
            [self._compare_distance(idx) for idx in range(2)]
        )
        df = pd.concat(dfs, axis=1)
        df.index = [self.target]  # type: ignore
        return df

##### Analysis for protein-ligand interactions (PLI) #####

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
        affinity_cols = ['score', 'sucos', 'centroid shift','shape_similarity',
                         'electrostatic_similarity', 'LE_heavyatom', 'LE_mw']
        boolean_cols = ['no_clashes', 'fully_matched', 'matched_rate']
        
        dfs = (
        [self._get_median_iqr(col) for col in affinity_cols] +
        [self._get_mean(col) for col in boolean_cols] +
        [self._count_interactions()] +
        [self._score_significance()]
        )
        df = pd.concat(dfs, axis=1)
        df.index = [self.target]  # type: ignore
        return df


class PLIAnalysisScore(PLIAnalysis):
    """Analysis for protein-ligand interactions (PLI) related metrics from scoring result of 3D *in-situ* methods.
    """
    def __init__(self, collections: DataCroupier):
        super().__init__(collections)
        self.test = collections.test_data.Score
        self.test_numerical, self.test_interaction = self._split_attr('test') # type: ignore

##### Analysis for property distributions #####

class PropAnalysis(AnalysisBase):
    """Analysis for property distributions."""
    def __init__(self, collections: DataCroupier):
        super().__init__(collections, 'Prop')
        self.props_croupier = [self._split_attr(attr) for attr in self.croupier_keys]
        self.test_desc, self.test_struc, self.alert = self.props_croupier[0] # type: ignore
        self.ref_desc = self.props_croupier[1][0] # type: ignore
        self.decoy_desc = self.props_croupier[2][0] # type: ignore

    def descriptor_dist(self) -> pd.DataFrame:
        
        keys = self.ref_desc.keys() # type: ignore
        results = []

        for key in keys:
            ref = self.ref_desc[key] # type: ignore
            test = self.test_desc[key] #type: ignore
            decoy = self.decoy_desc[key] #type: ignore

            w_ref, w_decoy, w_ref2decoy = wasserstein_ref_decoy(ref, decoy, test)
            w_shift = (w_ref + w_decoy) / w_ref2decoy

            ks_value, ks_significance = ks_distance(ref, test)

            results.append({
                'Descriptor': key,
                'Wasserstein_Ref': w_ref,
                'Wasserstein_Shift': w_shift,
                'KS_Ref': ks_value,
                'Significance': ks_significance <= 0.05
            })
        return pd.DataFrame(results)
    
    def analysis(self, descriptor_info:bool=True) -> tuple[tuple, pd.DataFrame]:
        """
        Returns:
            tuple[list[pd.Series], pd.DataFrame]: Series of mean values of descriptor distances, structural properties, and alert information,  
            and a DataFrame of descriptor distribution information.
        """
        desc_df = self.descriptor_dist()
        dfs = (
            desc_df.iloc[:, 1:].mean().to_frame().T,
            pd.DataFrame(self.test_struc).mean().to_frame().T,
            pd.DataFrame(self.alert).mean().to_frame().T
        )
        for df in dfs:
            df.index = [self.target] # type: ignore
        descri_info = desc_df if descriptor_info else pd.DataFrame()
        return dfs, descri_info
    

#TODO: Add cross-target analysis