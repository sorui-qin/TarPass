'''
Author: Rui Qin
Date: 2025-03-10 12:03:05
LastEditTime: 2025-05-14 21:41:59
Description: 
'''
import logging
import sys
from pathlib import Path
from tqdm import tqdm
from utils.constant import DASHLINE

class TqdmHandler(logging.StreamHandler):
    def emit(self, record):
        msg = self.format(record)
        tqdm.write(msg, end="\n")

def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.propagate = False # Disable propagation to avoid duplicate logs
    logger.setLevel(logging.DEBUG)

    # Console log
    console_handler = TqdmHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter(
        "%(levelname)s: %(message)s",
    )
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    return logger

project_logger = get_logger("project")

def log_config(project_logger, args):
    project_logger.info(DASHLINE)
    project_logger.info("\n" + "="*20 + " CONFIGURATION " + "="*20)
    for key, value in args.__dict__.items():
        project_logger.info(f"{key:10} : {str(value):<}")
    project_logger.info(DASHLINE)