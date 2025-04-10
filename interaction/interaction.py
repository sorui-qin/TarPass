'''
Author: Rui Qin
Date: 2025-04-10 20:57:37
LastEditTime: 2025-04-10 20:58:48
Description: 
'''
import re, json
from rdkit import Chem
from collections import defaultdict
from plip.structure.preparation import PDBComplex, PLInteraction
from utils.io import *
from utils.preprocess import read_in
from utils.logger import project_logger
from utils.constant import TARGETS, DASHLINE, INTERACTION_TYPES
from collections import defaultdict

def combined_pdb_complex(pdb:str, ligand:Chem.Mol, output_pdb):
    protein = Chem.MolFromPDBFile(pdb, sanitize=True)
    combined = Chem.CombineMols(protein, ligand)
    Chem.MolToPDBFile(combined, output_pdb)

def get_interactions(interactions: PLInteraction) -> defaultdict:
    check_list = [
        interactions.hbonds_ldon,
        interactions.hbonds_pdon,
        interactions.halogen_bonds,
        interactions.hydrophobic_contacts,
        interactions.metal_complexes,
        interactions.pication_laro,
        interactions.pication_paro,
        interactions.pistacking,
        interactions.saltbridge_lneg,
        interactions.saltbridge_pneg,
        #interactions.water_bridges, # not used
    ]

    interaction_types = INTERACTION_TYPES
    inters = defaultdict(set)
    for name, inter in zip(interaction_types, check_list):
        for i in inter:
            inters[name].add(f'{i.reschain}_{i.restype}{i.resnr}')
        inters[name] = sorted(inters[name], key=lambda x: int(re.search(r'\d+', x).group()))
    inters['H-Bond'] = inters['H-Bond acceptor'] + inters['H-Bond donor']
    return inters