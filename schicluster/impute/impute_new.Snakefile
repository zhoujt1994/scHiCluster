
# Load packages
import cooler
import pathlib
import pandas as pd

# Read and validate cell ids and get chromosome list

if 'cell_ids' not in locals():
    cell_table = pd.read_csv('cell_table.csv', index_col=0, header=None).squeeze(axis=1)
    cell_ids = cell_table.index

chromosomes = pd.read_csv(chrom_size_path, sep='\t', index_col=0, header=None).index

print(len(cell_ids), 'cells to process')
print(len(chromosomes), 'chromosomes in each cell.')

# Final targets
rule summary:
    input:
        expand('{cell_id}.cool', cell_id = cell_ids)
    shell:
        'touch Success && rm -rf impute_*_tmp && rm -rf .snakemake'

# Impute each chromosome of each cell
if 'input_scool' in locals():
    rule impute_chrom:
        input:
            input_scool
        output:
            temp('impute_{cell_id}_tmp/{chrom}.npz')
        threads:
            1
        shell:
            'hic-internal impute-chromosome '
            '--scool_url {input_scool}::/cells/{wildcards.cell_id} '
            '--chrom {wildcards.chrom} '
            '--resolution {resolution} '
            '--output_path {output} '
            '{logscale_str} '
            '--pad {pad} '
            '--std {std} '
            '--rp {rp} '
            '--tol {tol} '
            '--window_size {window_size} '
            '--step_size {step_size} '
            '--output_dist {output_dist} '
            '--min_cutoff {min_cutoff}'
elif 'cell_table' in locals():
    rule impute_chrom:
        output:
            temp('impute_{cell_id}_tmp/{chrom}.npz')
        params:
            contact_path=lambda wildcards: cell_table.loc[wildcards.cell_id]
        threads:
            1
        shell:
            'hic-internal impute-chromosome '
            '--contact_path {params.contact_path} '
            '--chrom_size_path {chrom_size_path} '
            '--chrom {wildcards.chrom} '
            '--resolution {resolution} '
            '--output_path {output} '
            '{logscale_str} '
            '--pad {pad} '
            '--std {std} '
            '--rp {rp} '
            '--tol {tol} '
            '--window_size {window_size} '
            '--step_size {step_size} '
            '--output_dist {output_dist} '
            '--min_cutoff {min_cutoff} '
            '--chr1 {chrom1} '
            '--chr2 {chrom2} '
            '--pos1 {pos1} '
            '--pos2 {pos2}'


# Aggregate chromosome HDF files for the same cells
rule agg_cell:
    input:
        expand('impute_{{cell_id}}_tmp/{chrom}.npz', chrom=chromosomes)
    output:
        '{cell_id}.cool'
    threads:
        1
    shell:
        'hic-internal aggregate-chromosomes '
        '--chrom_size_path {chrom_size_path} '
        '--resolution {resolution} '
        '--input_dir impute_{wildcards.cell_id}_tmp '
        '--output_path {output} '
        '--chrom_wildcard "{{chrom}}.npz"'
