'''
Author: Rui Qin
Date: 2025-06-26 19:02:36
LastEditTime: 2025-07-19 14:27:30
Description: 
'''
from typing import Literal
import numpy as np
import scikit_posthocs as sp
import statsmodels.stats.api as sms
from scipy import stats
from statsmodels.stats.multitest import multipletests

###### Statistical Tests #####

def test_normality(data) -> tuple[bool, float]:
    """Test if a sample follows a normal distribution using Shapiro-Wilk test or 
        D'Agostino and Pearson's test.
    """
    if len(data) < 20:
        _, p_value = stats.shapiro(data)
    else:
        _, p_value = stats.normaltest(data)
    is_normal = p_value > 0.05
    return is_normal, p_value


def test_variance_homogeneity(*data_groups) -> tuple[bool, float]:
    """Test if multiple samples have equal variances using Levene's test.
    Args:
        *data_groups: Variable number of data arrays to test.
    """
    _, p_value = stats.levene(*data_groups)
    is_equal_variance = p_value > 0.05
    return is_equal_variance, p_value

#### Significance Tests ####

# two groups
def ttest(
        data1, data2, equal_var:bool, 
        alternative:Literal['two-sided', 'less', 'greater'] = 'two-sided'
        ) -> tuple[bool, np.float64]:
    """Perform a t-test to compare the means of two independent samples.
    Args:
        data1 (array-like): First sample data.
        data2 (array-like): Second sample data.
        equal_var (bool): If True, perform a standard independent t-test; if False, perform Welch's t-test.
        alternative (str): Defines the alternative hypothesis. Options are 'two-sided', 'less', or 'greater'.
    """
    _, p_value = stats.ttest_ind(
        data1, data2, 
        equal_var=equal_var, 
        alternative=alternative
        )
    is_significant = p_value < 0.05 # type: ignore
    return is_significant, p_value # type: ignore

def mann_whitney_u(
        data1, data2, 
        alternative:Literal['two-sided', 'less', 'greater'] = 'two-sided'
        ) -> tuple[bool, float]:
    """Perform the Mann-Whitney U test to compare two independent samples.
    Args:
        data1 (array-like): First sample data.
        data2 (array-like): Second sample data.
        alternative (str): Defines the alternative hypothesis. Options are 'two-sided', 'less', or 'greater'.
    """
    _, p_value = stats.mannwhitneyu(
        data1, data2, 
        alternative=alternative
        )
    is_significant = p_value < 0.05
    return is_significant, p_value

# multiple groups
def anova(*data_groups, equal_var:bool) -> tuple[bool, float]:
    """Perform a one-way ANOVA test to compare means across multiple groups.
    Args:
        *data_groups: Variable number of data arrays to test.
        equal_var (bool): If True, assumes equal variances across groups; if False, uses Welch's ANOVA.
    """
    use_var = 'equal' if equal_var else 'unequal'
    p_value = sms.anova_oneway(data=data_groups, use_var=use_var).pvalue # type: ignore
    is_significant = p_value < 0.05
    return is_significant, p_value

def kruskal(*data_groups) -> tuple[bool, float]:
    """Perform the Kruskal-Wallis H-test for independent samples.
    Args:
        *data_groups: Variable number of data arrays to test.
    """
    _, p_value = stats.kruskal(*data_groups)
    is_significant = p_value < 0.05
    return is_significant, p_value

##### Post-hoc Tests #####

def dunn(
        *data_groups, control, 
        p_adjust:Literal['bonferroni', 'fdr_bh'] = 'fdr_bh'
        ) -> tuple[np.ndarray, np.ndarray]:
    """Perform Dunn's post-hoc test for multiple comparisons after Kruskal-Wallis test
    Args:
        *data_groups: Variable number of data arrays to test.
        control: The control group data array.
        p_adjust (str): Method for p-value adjustment. Options are 'bonferroni
            ' for Bonferroni correction or 'fdr_bh' for Benjamini-Hochberg FDR correction.
    """
    combine = [control] + [d for d in data_groups if d is not control]
    p_values= sp.posthoc_dunn(combine, p_adjust=p_adjust).values[0, 1:]
    is_significant = p_values < 0.05
    return is_significant, p_values

def dunnett(*data_groups, control) -> tuple[np.ndarray, np.ndarray]:
    """Perform Dunnett's test for multiple comparisons against a control group.
    Args:
        *data_groups: Variable number of data arrays to test.
        control: The control group data array.
        alternative (str): Defines the alternative hypothesis. Options are 'two-sided', 'less', or 'greater'.
    """
    p_values = stats.dunnett(*data_groups, control=control).pvalue
    is_significant = p_values < 0.05
    return is_significant, p_values

def tamhane(*data_groups, control) -> tuple[np.ndarray, np.ndarray]:
    """Perform Tamhane's T2 test for multiple comparisons.
    Args:
        *data_groups: Variable number of data arrays to test.
        control: The control group data array.
    """
    combine = [control] + [d for d in data_groups if d is not control]
    p_values = sp.posthoc_tamhane(combine, welch=True).values[0, 1:]
    is_significant = p_values < 0.05
    return is_significant, p_values


##### Multiple Testing Correction #####

def multiple_correction(p_values: np.ndarray, 
                        method:Literal['bonferroni', 'fdr_bh'] = 'fdr_bh'
                        ) -> tuple[np.ndarray, np.ndarray]:
    """Apply multiple testing correction to a set of p-values.
    Args:
        p_values (np.ndarray): Array of p-values to correct.
        method (str): Method for correction. Options are `bonferroni` for Bonferroni correction or \
            `fdr_bh` as default for Benjamini-Hochberg FDR correction.
    """
    rejected, corrected_pvalues, _, _ = multipletests(p_values, method=method)
    return rejected, corrected_pvalues

###### Effect Size #####

def cohen_d(data1, data2) -> tuple[float, str]:
    """Calculate Cohen's d effect size for two independent samples.  
    Interpretation is refered from Cohen, J. (1988). 
    Statistical power analysis for the behavioral sciences (2nd ed.). Hillside, NJ: Lawrence Erlbaum Associates.
    """
    nx, ny = len(data1), len(data2)
    s1, s2 = np.var(data1, ddof=1), np.var(data2, ddof=1)
    s_pooled = np.sqrt(((nx - 1)*s1 + (ny - 1)*s2) / (nx + ny - 2))
    d = (np.mean(data1) - np.mean(data2)) / s_pooled
    size = 'negligible' if abs(d) < 0.2 else 'small' if abs(d) < 0.5 else \
        'medium' if abs(d) < 0.8 else 'large'
    return d, size

def cliff_delta(data1, data2) -> tuple[float, str]:
    """Calculate Cliff's delta effect size for two independent samples.  
    Interpretation is refered from Romano, J., et al. (2006).
    Annual meeting of the Florida Association of Institutional Research. Vol. 177. No. 34
    """
    n1, n2 = len(data1), len(data2)
    diff = data1[:, None] - data2[None, :]
    count = np.sum(diff > 0) - np.sum(diff < 0)
    delta = count / (n1 * n2)
    size = 'negligible' if abs(delta) < 0.147 else 'small' if abs(delta) < 0.33 else \
        'medium' if abs(delta) < 0.474 else 'large'
    return delta, size

def omega_sq(*groups) -> tuple[float, str]:
    """Calculate the omega squared effect size for multiple groups.  
    Interpretation is refered from Field, A (2013) Discovering statistics using IBM SPSS Statistics. 
    Fourth Edition. Sage:London.
    """
    all_data = np.concatenate(groups)
    k = len(groups)
    n_total = len(all_data)
    grand_mean = np.mean(all_data)
    ss_between = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in groups)
    df_between = k - 1
    ss_within = sum(((g - np.mean(g))**2).sum() for g in groups)
    df_within = n_total - k
    ss_total = ((all_data - grand_mean)**2).sum()
    ms_within = ss_within / df_within

    omega_sq = (ss_between - df_between * ms_within) / (ss_total + ms_within)
    size = 'small' if omega_sq < 0.06 else 'medium' if omega_sq < 0.14 else 'large'
    return omega_sq, size

def epsilon_sq(*groups) -> tuple[float, str]:
    """Calculate the epsilon squared effect size for multiple groups.  
    Interpretation is refered from Rea, L. M., & Parker, R. A. (1992). Designing and conducting survey research:
    a comprehensive guide. San Francisco: Jossey-Bass Publishers.
    """
    h, _ = stats.kruskal(*groups)
    n = sum(len(g) for g in groups)
    k = len(groups)
    epsilon_sq = (h - k + 1) / (n - k)
    size = 'small' if epsilon_sq < 0.04 else 'medium' if epsilon_sq < 0.16 else 'large'
    return epsilon_sq, size