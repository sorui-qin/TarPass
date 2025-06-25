'''
Author: Rui Qin
Date: 2025-03-01 14:04:21
LastEditTime: 2025-06-25 15:14:31
Description: 
'''
import argparse
from .dockeval import dockeval_execute
from .moleeval import moleeval_execute
from .eval import eval_execute

def setup_arguments(parser: argparse.ArgumentParser):
    return parser

__all__ = ['setup_arguments', 'dockeval_execute', 'moleeval_execute', 'eval_execute']