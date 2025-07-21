'''
Author: Rui Qin
Date: 2025-07-07 17:22:34
LastEditTime: 2025-07-21 21:15:02
Description: 
'''
import ast
import numpy as np
import pandas as pd
from analysis_class import MolDisAnalysis, PLIAnalysis, PropAnalysis
from analysis.collect_eval import DataCroupier
from utils.constant import ROOT, TARGETS

DATA_DIR = ROOT / 'data'
REF_PKL = DATA_DIR / 'Reference_eval_results.pkl'
DEC_PKL = DATA_DIR / 'Decoy_eval_results.pkl'


def df_add_mean(df:pd.DataFrame) -> pd.DataFrame:
    """Add a row of mean values to the DataFrame."""
    df.index = TARGETS # type: ignore
    df.loc['Average'] = df.mean()
    return df

def triple_mean(series:pd.Series) -> str:
    triples = [ast.literal_eval(j) for j in series.values]
    return str([round(n,4) for n in np.mean(triples, axis=0)])

def process_sig(sigcol: pd.DataFrame) -> pd.DataFrame:
    def count_sig(df: pd.DataFrame) -> float:
        count = 0
        for _, series in df.iterrows():
            if series.iloc[2] and series.iloc[0] != 'negligible':
                count += 1
        return count / len(df)
    ref_sig = sigcol.iloc[:, :4]
    decoy_sig = sigcol.iloc[:, 4:]

    row = {}
    for i,col in enumerate(sigcol):
        row[col] = count_sig(ref_sig) if i == 3 \
            else count_sig(decoy_sig) if i == 7 \
            else 'NaN'
    sigcol.loc['Average'] = pd.Series(row).T
    return sigcol


def run_analysis(tests_data):
    results = {'pli': [], 'dist': [], 'stru': [], 'aler': []}
    for target in TARGETS:
        croupier = DataCroupier(tests[target], refs[target], decoys[target], target) # type: ignore

        prop_results = PropAnalysis(croupier).analysis()[0]
        mol_dist = MolDisAnalysis(croupier).analysis()
        pli_info = PLIAnalysis(croupier).analysis()

        desc_dist, stru_prop, alert_info = prop_results
        all_dist = pd.concat([mol_dist, desc_dist], axis=1)
        
        results['pli'].append(pli_info)
        results['dist'].append(all_dist)
        results['stru'].append(stru_prop)
        results['aler'].append(alert_info)

    results_concat = {k: pd.concat(v, ignore_index=True) for k, v in results.items()}