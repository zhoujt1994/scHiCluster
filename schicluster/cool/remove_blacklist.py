import pandas as pd
from collections import defaultdict
import pybedtools
from functools import lru_cache


@lru_cache()
def prepare_2d_blacklist_dict(blacklist_bedpe, resolution=10000):
    # read blacklist bed file, turn the region into bad pixel idx dict with resolution
    blacklist_bedpe_df = pd.read_csv(blacklist_bedpe, sep='\t', header=None)

    # turn region into region pixel idx
    blacklist_bedpe_df[1] //= resolution
    blacklist_bedpe_df[2] //= resolution
    # in case the region is smaller than resolution, add one so at least one pixel is bad
    blacklist_bedpe_df.loc[2, (blacklist_bedpe_df[2] -
                               blacklist_bedpe_df[1]) < 1] += 1
    blacklist_bedpe_df[4] //= resolution
    blacklist_bedpe_df[5] //= resolution
    # in case the region is smaller than resolution, add one so at least one pixel is bad
    blacklist_bedpe_df.loc[5, (blacklist_bedpe_df[5] -
                               blacklist_bedpe_df[4]) < 1] += 1

    chrom_pair_bad_points = defaultdict(set)
    for idx, row in blacklist_bedpe_df.iterrows():
        for i in range(row[1], row[2]):
            for j in range(row[4], row[5]):
                chrom_pair_bad_points[row[0], row[3]].add((i, j))

    # return a dict, key is chrom pair
    return chrom_pair_bad_points


def _is_2d_blacklist(row, blacklist_2d):
    chrom1, chrom2, pos1, pos2 = row
    judge = (pos1, pos2) in blacklist_2d[(chrom1, chrom2)]
    return judge


def filter_contacts(contact_path,
                    chrom_size_path=None,
                    blacklist_1d_path=None,
                    blacklist_2d_path=None,
                    remove_duplicates=True,
                    resolution_2d=10000,
                    chrom1=1,
                    pos1=2,
                    chrom2=5,
                    pos2=6):
    contacts = pd.read_csv(contact_path,
                           header=None,
                           sep='\t',
                           dtype={
                               chrom1: str,
                               pos1: int,
                               chrom2: str,
                               pos2: int
                           })
    print(f"{contact_path.split('/')[-1]}: {contacts.shape[0]} input contacts.")

    if remove_duplicates:
        # remove duplicates
        contacts.drop_duplicates(subset=[chrom1, pos1, chrom2, pos2], inplace=True)

    if blacklist_1d_path is not None:
        blacklist_bed_df = pd.read_csv(blacklist_1d_path, sep='\t', index_col=None, header=None)
        blacklist_bed = pybedtools.BedTool.from_dataframe(blacklist_bed_df).sort(g=chrom_size_path)

        # determine blacklist 1d (either side overlap with 1D blacklist)
        left_bed_df = contacts[[chrom1, pos1, pos1]].reset_index().iloc[:, [1, 2, 3, 0]]
        right_bed_df = contacts[[chrom2, pos2, pos2]].reset_index().iloc[:, [1, 2, 3, 0]]
        left_bed_df.columns = ['chrom', 'start', 'end', 'id']
        right_bed_df.columns = ['chrom', 'start', 'end', 'id']
        contact_bed_df = pd.concat([left_bed_df, right_bed_df])
        contact_bed = pybedtools.BedTool.from_dataframe(contact_bed_df).sort(g=chrom_size_path)
        # collect contact ids with either side overlap with blacklist
        bad_contacts = contact_bed.intersect(blacklist_bed, wa=True, u=True).to_dataframe()
        if bad_contacts.shape[0] > 0:
            bad_contacts = bad_contacts['name'].unique()
        else:
            # no bad contacts
            bad_contacts = set()
        # remove bad contacts
        contacts = contacts[~contacts.index.isin(bad_contacts)].copy()
        contact_bed.delete_temporary_history(ask=False)
        blacklist_bed.delete_temporary_history(ask=False)

    if blacklist_2d_path is not None:
        chrom_2d_blacklist = prepare_2d_blacklist_dict(blacklist_2d_path, resolution=resolution_2d)
        contacts_idx = contacts[[chrom1, chrom2, pos1, pos2]].copy()
        # turn contact location into bin idx with resolution
        contacts_idx[pos1] //= resolution_2d
        contacts_idx[pos2] //= resolution_2d

        # determine blacklist 2d (both side overlap with 2D blacklist)
        is_blacklist_2d = contacts_idx.apply(_is_2d_blacklist,
                                             blacklist_2d=chrom_2d_blacklist,
                                             axis=1)
        contacts = contacts[~is_blacklist_2d].copy()

    print(f"{contact_path.split('/')[-1]}: {contacts.shape[0]} filtered contacts in scool.")
    return contacts
