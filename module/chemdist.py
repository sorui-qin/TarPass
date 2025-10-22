'''
Author: Rui Qin
Date: 2025-07-19 14:27:02
LastEditTime: 2025-07-20 16:33:12
Description: 
'''
import numpy as np
from scipy import stats
from scipy.spatial.distance import jensenshannon
from sklearn.preprocessing import RobustScaler

def get_median_iqr(data):
    """Get median, IQR (Q1 & Q3) for the given data.
    Args:
        data (ArrayLike): Input data array.
    """
    data = np.array(data)[~np.isnan(data)]
    return (np.median(data),
            np.percentile(data, 25),
            np.percentile(data, 75))

def jsd(p, q) -> np.float64:
    """Calculate the Jensen-Shannon Divergence between two probability distributions.
    """
    return jensenshannon(p, q)

def ks_distance(p, q) -> tuple[np.float64, np.float64]:
    """Calculate the Kolmogorov-Smirnov distance and significance between two distributions.
    Ruturns:
        tuple: KS statistic and p-value.
    """
    ks_stat, p_value = stats.ks_2samp(p, q)
    return ks_stat, p_value # type: ignore

def wasserstein_norm(ref, test) -> np.float64:
    """Calculate the Wasserstein distance between reference and test distributions.
    """
    a_scale = np.array(ref).reshape(-1, 1)
    b_scale = np.array(test).reshape(-1, 1)
    scaler = RobustScaler()
    a_scale = scaler.fit_transform(a_scale).reshape(-1)
    b_scale = scaler.transform(b_scale).reshape(-1)
    return stats.wasserstein_distance(a_scale, b_scale)

def wasserstein_ref_rand(ref, rand, test) -> tuple[np.float64, np.float64, np.float64]:
    """Calculate the Wasserstein distance between test and reference or rand distributions.
    Args:
        ref (array-like): Reference distribution.
        rand (array-like): rand distribution.
        test (array-like): Test distribution.
    Returns:
        tuple: Wasserstein distances between reference and test, and rand and test distributions.
    """
    a_scale = np.array(ref).reshape(-1, 1)
    b_scale = np.array(rand).reshape(-1, 1)
    c_scale = np.array(test).reshape(-1, 1)
    scaler = RobustScaler()
    a_scale = scaler.fit_transform(a_scale).reshape(-1)
    b_scale = scaler.transform(b_scale).reshape(-1)
    c_scale = scaler.transform(c_scale).reshape(-1)
    return (
        stats.wasserstein_distance(a_scale, c_scale),
        stats.wasserstein_distance(b_scale, c_scale),
        stats.wasserstein_distance(a_scale, b_scale)
    )