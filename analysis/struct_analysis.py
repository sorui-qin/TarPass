'''
Author: Rui Qin
Date: 2025-07-10 15:07:11
LastEditTime: 2025-07-10 17:44:49
Description: 
'''
from typing import Literal, Optional
from collect_eval import MoleculesData, StructInfo, AnalysisBase
import numpy as np
import pandas as pd
from typing import cast
from module.significance import SignificanceTester


def get_median_iqr(data, prefix:Optional[str]=None) -> pd.DataFrame:
    """Get median, IQR (Q1 & Q3) for docking or scoring data."""
    prefix = f"{prefix}-" if prefix else ""
    stats = {
            f'{prefix}median': np.median(data),
            f'{prefix}q1': np.percentile(data, 25),
            f'{prefix}q3': np.percentile(data, 75)
        }
    return pd.DataFrame([stats])


