'''
Author: Rui Qin
Date: 2025-03-01 14:04:21
LastEditTime: 2025-06-14 17:19:10
Description: 
'''
import argparse
from .dockeval import dock_eval

def setup_arguments(parser: argparse.ArgumentParser):
    return parser

__all__ = ['setup_arguments', 'dock_eval']