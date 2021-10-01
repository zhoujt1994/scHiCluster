
# prepare variables
import pandas as pd
import pathlib
import cooler

cell_table = pd.read_csv(cell_table_path, index_col=0)
cell_ids = cell_table.index.tolist()
n_cells = len(cell_ids)
cell_id_to_url = cell_table['cell_url'].to_dict()

# get chromosomes
cool = cooler.Cooler(cell_table['cell_url'][0])
# instead of using cool.chromnames, using chrom_size_path chrom
chromnames = pd.read_csv(chrom_size_path, sep='\t', header=None, index_col=0, squeeze=True).to_dict()

# cleanup script
cleanup_cmd = 'rm -rf ' + ' '.join(chromnames)

matrix_types = ['E', 'E2', 'T', 'T2', 'Q']


# summary
rule summary:
    input:
        expand('{matrix_type}.cool', matrix_type=matrix_types)
    shell:
        cleanup_cmd + ' && touch finish'


# merge chromosomes into a single cool file
rule merge_chroms:
    input:
        expand('{chrom}.{{matrix_type}}.hdf', chrom=chromnames)
    output:
        '{matrix_type}.cool'
    threads:
        5
    shell:
        'hic-internal aggregate-chromosomes '
        '--chrom_size_path {chrom_size_path} '
        '--resolution {resolution} '
        '--input_dir ./ '
        '--output_path {output} '
        '--chrom_wildcard "{{chrom}}.{wildcards.matrix_type}.hdf"'


# merge cells' E matrix (npz) to group E and O matrix (coo stored in pd.HDFStore)
rule merge_E_O:
    input:
        expand('{chrom}/{cell_id}.E.npz', chrom=chromnames, cell_id=cell_ids)
    output:
        temp('{chrom}.E.hdf'),
        temp('{chrom}.E2.hdf')
    threads:
        5
    shell:
        'hic-internal merge-loop-matrix '
        '--output_dir {wildcards.chrom}/ '
        '--output_prefix {wildcards.chrom} '
        '--merge_type "E" '

# merge cells' T matrix (npz) to group T and T2 matrix (coo stored in pd.HDFStore)
rule merge_T_T2:
    input:
        expand('{chrom}/{cell_id}.T.npz', chrom=chromnames, cell_id=cell_ids)
    output:
        temp('{chrom}.T.hdf'),
        temp('{chrom}.T2.hdf')
    threads:
        5
    shell:
        'hic-internal merge-loop-matrix '
        '--output_dir {wildcards.chrom}/ '
        '--output_prefix {wildcards.chrom} '
        '--merge_type "T" '


# merge cells' Q matrix (imputed, before normalization, scool) to group Q matrix (coo stored in pd.HDFStore)
rule merge_Q:
    input:
        cell_table_path
    output:
        temp('{chrom}.Q.hdf')
    threads:
        5
    shell:
        'hic-internal merge-cell-impute-matrix '
        '--cell_urls_path {input} '
        '--chrom {wildcards.chrom} '
        '--output_prefix {wildcards.chrom}'

# Impute each chromosome of each cell
if keep_cell_matrix:
    rule loop_bkg_chrom:
        output:
            E='{chrom}/{cell_id}.E.npz',
            T='{chrom}/{cell_id}.T.npz'
        params:
            cell_url=lambda wildcards: cell_id_to_url[wildcards.cell_id],
            output_prefix='{chrom}/{cell_id}'
        threads:
            1
        shell:
            'hic-internal calculate-loop-matrix '
            '--cell_url {params.cell_url} '
            '--chrom {wildcards.chrom} '
            '--resolution {resolution} '
            '--output_prefix {params.output_prefix} '
            '--dist {dist} '
            '--cap {cap} '
            '--pad {pad} '
            '--gap {gap} '
            '--min_cutoff {min_cutoff} '
            '{log_e_str}'
else:
    rule loop_bkg_chrom:
        output:
            E=temp('{chrom}/{cell_id}.E.npz'),
            T=temp('{chrom}/{cell_id}.T.npz')
        params:
            cell_url=lambda wildcards: cell_id_to_url[wildcards.cell_id],
            output_prefix='{chrom}/{cell_id}'
        threads:
            1
        shell:
            'hic-internal calculate-loop-matrix '
            '--cell_url {params.cell_url} '
            '--chrom {wildcards.chrom} '
            '--resolution {resolution} '
            '--output_prefix {params.output_prefix} '
            '--dist {dist} '
            '--cap {cap} '
            '--pad {pad} '
            '--gap {gap} '
            '--min_cutoff {min_cutoff} '
            '{log_e_str}'
