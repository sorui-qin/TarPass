'''
Author: Rui Qin
Date: 2025-06-26 19:02:36
LastEditTime: 2025-06-27 12:43:43
Description: 
'''
import numpy as np
from scipy.spatial.distance import jensenshannon
from scipy.stats import wasserstein_distance, ks_2samp, entropy

def jsd(p, q):
    """Calculate the Jensen-Shannon Divergence between two probability distributions.
    """
    return jensenshannon(p, q)

def wasserstein(p, q):
    """Calculate the Wasserstein distance between two distributions.
    """
    return wasserstein_distance(p, q)