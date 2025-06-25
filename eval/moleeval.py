'''
Author: Rui Qin
Date: 2025-06-25 13:10:00
LastEditTime: 2025-06-25 15:12:27
Description: 
'''
import pandas as pd
from rdkit.Chem import Mol
from module.druglikeness import DruglikenessCalculator
from module.structural import StructuralCalculator
from utils.preprocess import read_in, to_mols
from utils.constant import Path, DASHLINE, TARGETS
from tqdm.contrib.concurrent import process_map
from utils.logger import project_logger


class MoleEval:
    """Evaluate molecules for druglikeness and structural properties.
    
    Args:
        target_path (str): Path to the target directory containing molecule files.
    """
    def __init__(self, smis: list[str]):
        self.smis = smis
        self.mols = to_mols(self.smis)

    def evaluate(self) -> list[dict]:
        """Evaluate all molecules in the target directory."""
        props = process_map(
            all_properties, self.mols, chunksize=100,
            desc=f"Evaluating molecular properties",
            )
        
        updated_props = []
        for i, smi in enumerate(self.smis):
            info = {'index': i, 'smi': smi}
            info.update(props[i])
            updated_props.append(info)
        return updated_props


def all_properties(mol: Mol) -> dict:
    """
    Calculate all properties of a molecule, including druglikeness and structural properties.
    Args:
        mol (Chem.rdchem.Mol): Single molecule object from RDKit.
    
    Returns:
        dict: A dictionary containing various properties of the molecule.
    """
    druglikeness = DruglikenessCalculator().calculate_all(mol)
    structural = StructuralCalculator().calculate_all(mol)
    properties = {**druglikeness, **structural}
    return properties


def mole_eval(smis: list[str]) -> list[dict]:
    """
    Evaluate a list of SMILES strings for druglikeness and structural properties.
    
    Args:
        smis (list[str]): List of SMILES strings to evaluate.
    
    Returns:
        list[dict]: List of dictionaries containing properties for each molecule.
    """
    evaluator = MoleEval(smis)
    return evaluator.evaluate()

############## Execution Functions ##############

def moleeval_execute(args):
    work_dir = Path(args.path)
    for target in TARGETS:

        # Check if target directory exists
        project_logger.info(DASHLINE)
        target_dir = work_dir / target
        if not target_dir.exists():
            project_logger.warning(f"Target folder {target_dir} not found, skipping evaluation.")
            continue

        # Check if results exists
        project_logger.info(f'Start evaluating molecules for {target}...')
        results_dir = target_dir / 'results'
        eval_output = results_dir / f'mole_eval_results.csv'
        if eval_output.exists():
            project_logger.info(f"Evaluation results already exist for {target}, skipping evaluation.")
            continue

        # Execute evaluation and save results
        project_logger.info(DASHLINE)
        smis, _ = read_in(target_dir)
        result = mole_eval(smis)
        pd.DataFrame(result).to_csv(eval_output, index=False)
        project_logger.info(f"Evaluation results saved to {eval_output}.")
        project_logger.info(DASHLINE)