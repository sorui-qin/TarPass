'''
Author: Rui Qin
Date: 2025-07-07 17:22:34
LastEditTime: 2025-07-22 19:25:39
Description: 
'''
import argparse
import ast
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
from analysis.analysis_class import (MolDisAnalysis, PLIAnalysis,
                                     PLIAnalysisScore, PropAnalysis)
from analysis.collect_eval import DataCroupier, MoleculesData, collect_eval
from utils.constant import ROOT, TARGETS
from utils.io import read_pkl, write_pkl
from utils.logger import log_config, project_logger

DATA_DIR = ROOT / 'data'
REF_PKL = DATA_DIR / 'Reference_eval_results.pkl'
DEC_PKL = DATA_DIR / 'Decoy_eval_results.pkl'

def df_add_mean(df:pd.DataFrame) -> pd.DataFrame:
    """Add a row of mean values to the DataFrame."""
    df.loc['Average'] = df.mean()
    return df

def triple_mean(series:pd.Series) -> str:
    """Calculate the mean of a series of triples and return as a formatted string."""
    triples = [ast.literal_eval(j) for j in series.values]
    return str([round(n,4) for n in np.mean(triples, axis=0)])

def count_sig(df:pd.DataFrame, effect_col=1, pval_col=3) -> float:
    return (
        (df.iloc[:, pval_col] & (df.iloc[:, effect_col] != 'negligible'))
        .sum() / len(df)
    )

def process_sig(sigcol:pd.DataFrame) -> pd.DataFrame:
    assert sigcol.shape[1] == 8, "sigcol should have 8 columns"
    sigcol = sigcol.copy()
    ref_sig = sigcol.iloc[:, :4]
    decoy_sig = sigcol.iloc[:, 4:]
    
    avg_row = {}
    for i, col in enumerate(sigcol.columns):
        if i == 3:  # ref significance column
            avg_row[col] = count_sig(ref_sig)
        elif i == 7:  # decoy significance column
            avg_row[col] = count_sig(decoy_sig)
        else:
            avg_row[col] = 'NaN'
    
    sigcol.loc['Average'] = pd.Series(avg_row)
    return sigcol

def process_pli(df_pli: pd.DataFrame) -> pd.DataFrame:
    """Process the PLI DataFrame to calculate means and significance."""
    assert df_pli.shape[1] == 19, "df_pli should have 19 columns"
    df_pli = df_pli.copy()
    
    tricol = df_pli.iloc[:, 0:7].copy()
    meancol = df_pli.iloc[:, 7:11].copy()
    sigcol = df_pli.iloc[:, 11:19].copy()
    
    tricol.loc['Average'] = tricol.apply(triple_mean)
    meancol = df_add_mean(meancol)
    sigcol = process_sig(sigcol)

    return pd.concat([tricol, meancol, sigcol], axis=1)

def _load_test(path: Path):
    if path.is_file() and path.suffix == '.pkl':
        return read_pkl(path), 'pkl'
    if path.is_dir():
        return {}, 'dir'
    raise ValueError(f"Invalid path: {path}. It should be directory or a pickle file.")

############## Execution Functions ##############

def setup_arguments(parser: argparse.ArgumentParser):
    parser.add_argument('--output', '-o', type=str, help='Directory of the output file where the results will be saved. \
                        Default same to the working directory.')
    return parser

def execute(args):
    log_config(project_logger, args)
    work_dir = Path(args.path)
    tests_data, read_source = _load_test(work_dir)
    
    refs = read_pkl(REF_PKL)
    decoys = read_pkl(DEC_PKL)

    results = {'pli_dock': [], 'pli_score': [],'dist': [], 'stru': [], 'aler': []}
    external_info = {'smis': {}, 'desc_details': {}}
    
    for target in tqdm(TARGETS, desc='Analysis targets'):
        # Extract the test data for the target
        if read_source == 'pkl':
            test_target = tests_data.get(target, None) # type: ignore
            if test_target is None:
                continue
        elif read_source == 'dir':
            test_target = collect_eval(work_dir / target / 'results')
        assert isinstance(test_target, MoleculesData), f"Invalid data for target {target}"

        croupier = DataCroupier(test_target, refs[target], decoys[target], target) # type: ignore

        prop_results, desc_detail = PropAnalysis(croupier).analysis()
        if croupier.test_data.Score:
            pli_score = PLIAnalysisScore(croupier).analysis()
            results['pli_score'].append(pli_score)
        else:
            results['pli_score'].append(pd.DataFrame())
        mol_dist = MolDisAnalysis(croupier).analysis()
        pli_info = PLIAnalysis(croupier).analysis()

        external_info['smis'][target] = croupier.test_data.Mol.smiles
        external_info['desc_details'][target] = desc_detail

        desc_dist, stru_prop, alert_info = prop_results
        all_dist = pd.concat([mol_dist, desc_dist], axis=1)
        
        results['pli_dock'].append(pli_info)
        results['dist'].append(all_dist)
        results['stru'].append(stru_prop)
        results['aler'].append(alert_info)

    results_concat = {k: pd.concat(v, ignore_index=False) for k, v in results.items()}

    # PLI Analysis
    results_concat['pli_dock'] = process_pli(results_concat['pli_dock'])
    if not results_concat['pli_score'].empty:
        results_concat['pli_score'] = process_pli(results_concat['pli_score'])

    # Other row
    for k in ['dist', 'stru', 'aler']:
        results_concat[k] = df_add_mean(results_concat[k])

    # TODO: cross-target analysis

    output_dir = Path(args.output) if args.output else work_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    for k, v in results_concat.items():
        if not v.empty:
            v.to_csv(output_dir / f'tarpass_anlysis_{k}.csv', index=True)
    write_pkl(output_dir / 'descriptor_details.pkl', external_info['desc_details'])
    project_logger.info(f"Analysis results were saved to {output_dir}.")