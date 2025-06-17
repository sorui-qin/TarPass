'''
Author: Rui Qin
Date: 2025-03-01 14:04:21
LastEditTime: 2025-06-17 17:18:53
Description: 
'''
import argparse
from .dockeval import dockeval_execute

def setup_arguments(parser: argparse.ArgumentParser):
    return parser

__all__ = ['setup_arguments', 'dockeval_execute']