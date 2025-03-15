'''
Author: Rui Qin
Date: 2025-03-15 13:52:13
LastEditTime: 2025-03-15 16:21:53
Description: 
'''
from utils.io import *
from utils.dock import LigPrep
from utils.logger import project_logger
from dock.docking_vina import VinaDock
from dock.docking_gnina import GninaDock


class Dock():
    def  __init__(self, mols, target, method='gnina', mode='dock'):
        self.mols = mols
        self.target = target
        self.mode = mode
        self.method = VinaDock if method == 'vina' else GninaDock
        pass


def setup_arguments(parser):
    parser.add_argument()


def execute(args):
    
    # Ligand preparation

    # Docking with AutoDock-Vina

    # Docking with Gnina