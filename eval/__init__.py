'''
Author: Rui Qin
Date: 2025-03-01 14:04:21
LastEditTime: 2025-06-18 16:34:54
Description: 
'''
import argparse
from .dockeval import dockeval_execute

def setup_arguments(parser: argparse.ArgumentParser):
    parser.add_argument('--format', default='json', type=str, choices=['csv', 'json'], 
                    help='Output format for the evaluation results. Default in json.')
    return parser

__all__ = ['setup_arguments', 'dockeval_execute']