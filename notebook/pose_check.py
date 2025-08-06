'''
Author: Rui Qin
Date: 2025-08-06 11:44:46
LastEditTime: 2025-08-06 16:40:37
Description: 
'''
import sys
import pandas as pd
from pathlib import Path
from posebusters import PoseBusters
from tqdm import tqdm
from utils.constant import TARGETS
from utils.preprocess import read_in
from utils.logger import disable_logger

disable_logger()

path = Path(sys.argv[1])
pb = PoseBusters('mol')

results = []
for target in tqdm(TARGETS, desc='Checking conformations'):
    target_dir = path / target
    if not target_dir.exists():
        print(f"{target} directory does not exist. Skipping...")
        continue
    _, mols = read_in(target_dir)
    df = pb.bust(mol_pred=mols, full_report=True)
    df.index = target # type: ignore
    results.append(df)

pd.concat(results, ignore_index=False).to_csv(path / 'pose_check.csv', index=True)