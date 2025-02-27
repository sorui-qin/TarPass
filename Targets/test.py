'''
Author: Rui Qin
Date: 2025-02-22 22:14:22
LastEditTime: 2025-02-24 16:57:57
Description: 
'''
from plip.structure.preparation import PDBComplex
my_mol = PDBComplex()
my_mol.load_pdb('7s1s_obabel.pdb')
my_bsid = '85K:A:702'
my_mol.analyze()
my_interactions = my_mol.interaction_sets[my_bsid] 
my_interactions