import numpy as np
import pandas as pd
from scipy.stats import f as f_dist


def one_way_anova(chrom_loop_ds, da_name, value_type, group_n_dim='group_n', group_dim='sample_id'):
    """
    Perform one-way ANOVA on a single-chrom loop dataset.

    Parameters
    ----------
    chrom_loop_ds
        A single-chrom loop dataset.
    da_name
        The name of the data array to perform ANOVA on.
    value_type
        The value type of the data array to perform ANOVA on.
        both "{value_type}" and "{value_type}2" should be present in the "{da_name}_value_type" dimension.
    group_n_dim
        The name of the group number variable.
    group_dim
        The name of the group dimension.

    Returns
    -------
    F statistics and P-values of the one-way ANOVA.
    """
    loop_ds = chrom_loop_ds
    # number of cells per group
    n = loop_ds[group_n_dim]
    n_group = n.size
    total_n = int(n.sum())

    # sum(x) / ni, ni is number of cells in i group
    x = loop_ds[da_name].sel({f'{da_name}_value_type': value_type})
    # sum(x^2) / ni
    x2 = loop_ds[da_name].sel({f'{da_name}_value_type': f'{value_type}2'})

    # total variance
    x_sum = (x * n).sum(dim=group_dim)
    x2_sum = (x2 * n).sum(dim=group_dim)
    sst = x2_sum - np.power(x_sum, 2) / total_n

    # within group variance
    ssw = ((x2 - np.power(x, 2)) * n).sum(dim=group_dim)

    # ssb between group variance; ssb = sst - ssw
    # ssb / ssw = sst / ssw - 1
    # f = (sst / ssw - 1) * (N - k) / (k - 1)
    f = (sst / ssw - 1) * (total_n - n_group) / (n_group - 1)
    f = f.to_pandas()

    # p value from F distribution
    p = f_dist(n_group - 1, total_n - n_group).sf(f)
    p = pd.Series(p, index=f.index)
    return f, p


def merge_groups(loop_ds, group_map, da_name, group_dim='sample_id', group_n_dim='group_n'):
    """
    Merge groups into larger groups in a loop dataset.

    Parameters
    ----------
    loop_ds
        A loop dataset.
    group_map
        A pd.Series mapping from old group names to new group names.
    da_name
        The name of the data array to merge groups for.
    group_dim
        The name of the group dimension.
    group_n_dim
        The name of the group number variable.

    Returns
    -------
    A loop dataset with merged groups.
    """
    loop_ds[da_name] = loop_ds[da_name] * loop_ds.coords[group_n_dim]

    loop_ds['_sample_group'] = group_map
    loop_ds = loop_ds.groupby('_sample_group').sum(dim=group_dim)

    cell_count = loop_ds.coords[group_n_dim].to_pandas()
    sample_group_count = cell_count.groupby(group_map).sum()
    sample_group_count.index.name = '_sample_group'
    loop_ds.coords[group_n_dim] = sample_group_count

    loop_ds[da_name] = loop_ds[da_name] / loop_ds[group_n_dim]

    loop_ds = loop_ds.rename({
        '_sample_group': group_dim
    })
    return loop_ds
