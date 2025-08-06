'''
Author: Rui Qin
Date: 2025-08-06 11:44:46
LastEditTime: 2025-08-06 20:50:55
Description: 
'''
import sys
import pandas as pd
from pathlib import Path
from posebusters import PoseBusters
from tqdm import tqdm
from utils.constant import TARGETS
from utils.io import write_pkl
from utils.preprocess import read_in
from utils.logger import disable_logger

disable_logger()

path = Path(sys.argv[1])
pb = PoseBusters('mol')

results = []
store = {}
for target in tqdm(TARGETS, desc='Checking conformations'):
    target_dir = path / target
    if not target_dir.exists():
        print(f"{target} directory does not exist. Skipping...")
        continue
    _, mols = read_in(target_dir)
    df = pb.bust(mol_pred=mols, full_report=True)
    store[target] = df
    series = df.iloc[:, :11].mean(skipna=True)
    series.name = target
    results.append(series)

all_df = pd.concat(results, axis=1).T
all_df.loc['Average'] = all_df.mean()
all_df.to_csv(path / 'pose_check.csv', index=True)
write_pkl(path / 'pose_check_details.pkl', store)