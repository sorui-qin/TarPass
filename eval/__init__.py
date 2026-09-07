"""
Author: Rui Qin
Date: 2025-03-01 14:04:21
LastEditTime: 2025-07-07 16:12:22
Description:
"""

import argparse

from .dockeval import dockeval_execute
from .eval_main import eval_execute
from .moleeval import moleeval_execute


def setup_arguments(parser: argparse.ArgumentParser):
    return parser


__all__ = ["dockeval_execute", "eval_execute", "moleeval_execute", "setup_arguments"]
