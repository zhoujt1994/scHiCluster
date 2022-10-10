import numpy as np
from scipy.stats import chi2_contingency


def diff_bound(boundary_count, group_n):
    no_boundary = group_n[:, None] - boundary_count
    stats = np.zeros(boundary_count.shape[1])
    pv = np.ones(boundary_count.shape[1])
    bin_filter = np.logical_and(
        boundary_count.sum(axis=0) > 0,
        no_boundary.sum(axis=0) > 0
    )

    for i in range(boundary_count.shape[1]):
        if bin_filter[i]:
            contingency = [boundary_count[:, i], no_boundary[:, i]]
            stats[i], pv[i], _, _ = chi2_contingency(contingency)
    return stats, pv
