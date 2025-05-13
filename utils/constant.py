'''
Author: Rui Qin
Date: 2025-03-15 16:12:18
LastEditTime: 2025-05-13 17:32:24
Description: 
'''
from pathlib import Path

ROOT = Path(__file__).parent.parent
TMP = ROOT / 'tmp'

TARGETS = [
    '5HT2A',
    '5HT2A-AF',
    'BCL2',
    'BRAF',
    'BRD4',
    'BRD4-holo',
    'BTK',
    'Beta2AR',
    'DRD2',
    'HDAC6',
    'HIV-RT',
    'JAK2',
    'MEK1',
    'NAMPT',
    'PI3K-alpha',
    'PPAR-alpha',
    'PRMT5',
    'ROCK1',
    'RXR-alpha',
    'TYK2'
    ]

DASHLINE = '='*50

INTERACTION_TYPES = [
    'H-Bond acceptor',
    'H-Bond donor',
    'Halogen bond',
    'Hydrophobic contact',
    #'Metal complexes',
    'Pi-cation lAro',
    'Pi-cation pAro',
    'Pi-stacking',
    'Salt Bridge lneg',
    'Salt Bridge pneg',
    #'Water bridge'
    ]