'''
Author: Rui Qin
Date: 2025-05-13 23:35:54
LastEditTime: 2025-05-14 16:46:30
Description: 
'''
import argparse
import json
from interaction_tools import *
from multiprocessing import Pool
from tqdm import tqdm
from functools import partial
from utils.io import temp_manager, read_sdf
from utils.constant import ROOT, TARGETS

def analyze_tmppdb(mol:Chem.Mol, empty_pdb:Path) -> defaultdict:
    """
    Analyze the interactions of a ligand with a target protein using PLIP.
    Args:
        mol (Chem.Mol): The ligand molecule.
        empty_pdb (path): Path to the empty PDB file of corresponding target.
    Returns:
        defaultdict: Interactions categorized by type.
    """
    with temp_manager('pdb') as tmpdir:
        combined_pdb_complex(empty_pdb, mol, tmpdir)
        interactions = analyze_plip(tmpdir)
        return interactions

def setup_arguments(parser: argparse.ArgumentParser):
    parser.add_argument('-p', '--path', required=True, type=str, help='path to the folder where generated molecules for testing will be stored.')
    parser.add_argument('-m', '--mode', required=True, type=str, choices=['dock', 'generate'], help='Anlysis interactions with docking pose or generated pose.')
    return parser

def execute(args):
    plip_tmp()
    for target in TARGETS:
        empty_pdb = next((ROOT / f'Targets/{target}').glob(f'{target}*.pdb'))
        analyze_func = partial(analyze_tmppdb, empty_pdb=empty_pdb)
        with Pool() as pool:
            interactions = list(
                tqdm(pool.imap(analyze_func, mols),
                    desc=f"Analyzing interactions with {target}", total=len(mols)))
        

if __name__ == '__main__':
    args = setup_arguments(argparse.ArgumentParser()).parse_args()

    # Preparation for reading poses
    mode = args.mode
    pattern = ['result/gnina-dock_docking_results.pkl',
               'result/vina-dock_docking_results.pkl']

    # Analyze
    plip_tmp()
    allkey_inters = json.load(open(ROOT/'interaction/key_interactions.json'))

    for target in TARGETS:
        key_inters = allkey_inters[target]
        empty_pdb = next((ROOT / f'Targets/{target}').glob(f'{target}*.pdb'))

        mols = []
        # Running interaction analysis in parallel
        analyze_func = partial(analyze_tmppdb, empty_pdb=empty_pdb)
        with Pool() as pool:
            results = list(tqdm(pool.imap(analyze_func, mols),
                                desc="Analyzing interactions", total=len(mols)))
        