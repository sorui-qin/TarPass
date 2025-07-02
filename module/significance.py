'''
Author: Rui Qin
Date: 2025-07-02 16:59:50
LastEditTime: 2025-07-02 19:52:45
Description: 
'''
import numpy as np
import pandas as pd
from typing import Union
from utils.stats import (
    test_normality, test_variance_homogeneity,
    ttest, mann_whitney_u, anova, kruskal,
    dunn, dunnett, tamhane,
    cohen_d, cliff_delta, omega_sq, epsilon_sq
)

def washing_data(*data_groups) -> list[np.ndarray]:
    """Clean and validate input data groups.
    """
    cleaned = []
    for i, data in enumerate(data_groups):
        data_array = np.asarray(data)

        # Remove NaN and Inf values
        valid_mask = ~(np.isnan(data_array) | np.isinf(data_array))
        data_array = data_array[valid_mask]

        # Execute validation checks
        if data_array.size == 0:
            raise ValueError(f"Data group {i+1} is empty after cleaning.")
        if data_array.ndim > 1:
            raise ValueError(f"Data group {i+1} should be 1D, got {data_array.ndim}D.")
        if data_array.dtype.kind not in 'biufc':
            raise ValueError(f"Data group {i+1} must be numeric, got {data_array.dtype}.")
        if data_array.size < 3:
            raise ValueError(f"Data group {i+1} has insufficient sample size: {data_array.size}")

        cleaned.append(data_array)
    return cleaned


class SignificanceTester:
    """Statistical significance tester with automatic method selection.
    Args:
        control: Control group data
        data_groups: Test groups (list of arrays or single array)
    """
    
    def __init__(self, control, data_groups: Union[list, np.ndarray], metrics_name: str = 'default'):
        # Clean control group
        self.control = washing_data(control)[0]
        # Handle different input formats for data_groups
        if not isinstance(data_groups, (list, tuple)):
            data_groups = [data_groups]
        self.data_groups = washing_data(*data_groups)
        self.n_groups = len(self.data_groups)
        self.all_data = [self.control] + self.data_groups
        self.metrics_name = metrics_name

    def normality(self) -> tuple[bool, list[float]]:
        """Test normality for all groups.
        
        Returns:
            Tuple of (all_normal, detailed_results)
        """
        p, b = zip(*[test_normality(data) for data in self.all_data])
        return all(p), list(b)
    
    def variance_homogeneity(self) -> tuple[bool, float]:
        """Test if data groups have equal variance.
        """
        return test_variance_homogeneity(*self.all_data)
    
    def compare_two_groups(self, group_index: int = 0) -> dict:
        """
        Compare control with a specific test group.
        
        Args:
            group_index: Index of test group to compare (0-based)
        """
        if group_index >= self.n_groups:
            raise IndexError(f"Group index {group_index} out of range (0-{self.n_groups-1})")
        
        test_group = self.data_groups[group_index]
        is_normal, normal_p = self.normality()

        if is_normal:
            # Use parametric test
            equal_var, _ = self.variance_homogeneity()
            is_sig, p_val = ttest(self.control, test_group, equal_var=equal_var)
            test_name = f"{'Student' if equal_var else 'Welch'} t-test"
            effect_size, effect_interp = cohen_d(self.control, test_group)
        else:
            # Use non-parametric test
            is_sig, p_val = mann_whitney_u(self.control, test_group)
            test_name = "Mann-Whitney U test"
            effect_size, effect_interp = cliff_delta(self.control, test_group)
        
        return {
            'metric_name': self.metrics_name,
            'control_normal': normal_p[0],
            'test_normal': normal_p[group_index + 1],
            'equal_variance': locals().get('equal_var', None),
            'test_name': test_name,
            'p_value': p_val,
            'significant': is_sig,
            'effect_size': effect_size,
            'effect_interpretation': effect_interp,
        }
    

    def compare_multiple_groups(self) -> list[dict]:
        """
        Perform overall test for multiple groups.
        """
        if self.n_groups < 2:
            raise ValueError("Need at least 2 test groups for multiple comparison")
        
        # Check assumptions
        all_normal, normal_p = self.normality()
        
        # Choose appropriate test
        if all_normal:
            # Use ANOVA
            equal_var, _ = self.variance_homogeneity()
            is_sig, p_val = anova(*self.all_data, equal_var=equal_var)
            test_name = f"{'One-way' if equal_var else 'Welch'} ANOVA"
            effect_size = omega_sq(*self.all_data)
            if p_val < 0.05:
                if equal_var:
                    posthoc_name = "Dunnett's test"
                    posthoc_sigs, posthoc_ps = dunnett(*self.data_groups, control=self.control)
                else:
                    posthoc_name = "Tamhane's T2"
                    posthoc_sigs, posthoc_ps = tamhane(*self.data_groups, control=self.control)
                posthoc_effect_size, posthoc_effect_interp = zip(*[
                    cohen_d(self.control, group) for group in self.data_groups
                ])
        else:
            # Use Kruskal-Wallis
            is_sig, p_val = kruskal(*self.all_data)
            test_name = "Kruskal-Wallis H-test"
            effect_size = epsilon_sq(*self.all_data)
            if p_val < 0.05:
                posthoc_name = "Dunn's test"
                posthoc_sigs, posthoc_ps = dunn(*self.data_groups, control=self.control)
                posthoc_effect_size, posthoc_effect_interp = zip(*[
                    cliff_delta(self.control, group) for group in self.data_groups
                ])
        
        reports = []
        for i in range(self.n_groups):

            part1 = {
                'metric_name': self.metrics_name,
                'index': i,
                'total_groups': self.n_groups,
                'all_normal': all_normal,
                'normal_p_values': normal_p[i+1],
                'equal_variance': locals().get('equal_var', None),
                'test_name': test_name,
                'p_value': p_val,
                'significant': is_sig,
                'effect_size': effect_size,
            }

            posthoc = {} if p_val >= 0.05 else {
                'posthoc_name': posthoc_name,
                'posthoc_significant': posthoc_sigs[i],
                'posthoc_p_value': posthoc_ps[i],
                'posthoc_effect_size': posthoc_effect_size[i],
                'posthoc_effect_interpretation': posthoc_effect_interp[i]
            }
            
            reports.append({**part1, **posthoc})
        return reports
    
    def auto_analysis(self) -> list[dict]:
        """Automatically select and perform appropriate tests based on data characteristics.
        """
        if self.n_groups == 1:
            return [self.compare_two_groups(0)]
        else:
            return self.compare_multiple_groups()