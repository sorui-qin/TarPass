'''
Author: Rui Qin
Date: 2025-06-25 15:13:31
LastEditTime: 2025-06-25 17:10:14
Description: 
'''
import pandas as pd
from eval.dockeval import dock_eval
from eval.moleeval import mole_eval
from utils.constant import DASHLINE, TARGETS, Path
from utils.io import dump_json
from utils.logger import log_config, project_logger
from utils.preprocess import read_in


def eval_execute(args):
    log_config(project_logger, args)
    work_dir = Path(args.path)
    for target in TARGETS:

        # Check if target directory exists
        target_dir = work_dir / target
        if not target_dir.exists():
            project_logger.warning(f"Target folder {target_dir} not found, skipping evaluation.")
            continue
        
        # Reading...
        results_dir = target_dir / 'results'
        project_logger.info(DASHLINE)
        smis, mols = read_in(target_dir)

        # Docking evaluation
        dock_output = results_dir / f'dock_eval_results.json'
        if not dock_output.exists():
            project_logger.info(f'Start evaluating docking results for {target}...')
            dock_res = dock_eval(mols, target_dir)
            dump_json(dock_output, dock_res)
            project_logger.info(f"Evaluation results saved to {dock_output}.")
        else:
            project_logger.info(f"Docking valuation results already exist for {target}, skipping to the next...")

        # Molecule evaluation
        mole_output = results_dir / f'mole_eval_results.csv'
        if not mole_output.exists():
            project_logger.info(f'Start evaluating molecules for {target}...')
            mole_res = mole_eval(smis)
            pd.DataFrame(mole_res).to_csv(mole_output, index=False)
        else:
            project_logger.info(f"Evaluation results already exist for {target}, skipping evaluation.")

        project_logger.info(DASHLINE)