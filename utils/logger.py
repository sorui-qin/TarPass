'''
Author: Rui Qin
Date: 2025-03-10 12:03:05
LastEditTime: 2025-03-10 14:08:23
Description: 
'''
import logging
import sys
from pathlib import Path
from tqdm import tqdm

class TqdmHandler(logging.StreamHandler):
    def emit(self, record):
        msg = self.format(record)
        tqdm.write(msg, end="\n")

def get_logger(name: str, log_file: str = None) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # Console log
    console_handler = TqdmHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter(
        "%(levelname)s: %(message)s",
    )
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # File log
    #if log_file:
        #log_path = Path(log_file).parent
        #log_path.mkdir(parents=True, exist_ok=True)
        
        #file_handler = logging.FileHandler(log_file)
        #file_handler.setLevel(logging.DEBUG)
        #file_formatter = logging.Formatter(
            #"%(levelname)s - %(message)s")
        #file_handler.setFormatter(file_formatter)
        #logger.addHandler(file_handler)

    return logger

project_logger = get_logger("project")