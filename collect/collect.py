'''
Author: Rui Qin
Date: 2025-07-11 11:42:13
LastEditTime: 2025-07-11 12:06:16
Description: 
'''
import argparse
from pathlib import Path
from analysis.collect_eval import collect_eval_all
from utils.logger import log_config, project_logger


def setup_arguments(parser: argparse.ArgumentParser):
    parser.add_argument('--output', '-o', type=str, help='Directory of the output file where the results will be saved. \
                        Default same to the working directory.')
    parser.add_argument('--prefix', '--pre', type=str, help='Prefix used for naming the output file. \
                        Default same to the name of the working directory.')
    return parser


def execute(args):
    log_config(project_logger, args)
    path = Path(args.path)
    if args.prefix is None:
        args.prefix = path.name
    collect_eval_all(path, args.output, args.prefix)