'''
Author: Rui Qin
Date: 2025-03-17 20:15:43
LastEditTime: 2025-03-17 20:18:36
Description: 
'''
from utils.io import *
from utils.preprocess import read_in
from utils.logger import project_logger
from utils.constant import TARGETS, DASHLINE
from rdkit import Chem

def combined_pdb_complex(pdb:str, ligand:Chem.Mol, output_pdb):
    protein = Chem.MolFromPDBFile(pdb)
    combined = Chem.CombineMols(protein, ligand)
    Chem.MolToPDBFile(combined, output_pdb)