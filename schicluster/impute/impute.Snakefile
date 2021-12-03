
# config example, this will be added by code
# input_scool = '../contacts/full.10K.scool'
# chrom_size_path = '/home/hanliu/ref/mouse/genome/mm10.main.nochrM.chrom.sizes'
# logscale = False
# pad = 1
# std = 1
# window_size = int(4e7)
# step_size = int(1e7)
# output_dist = int(1e7)
# resolution = 10000
# rp = 0.5
# tol = 0.01
# min_cutoff = 1e-5
# cell_ids = ['CEMBA191126_9J_1-CEMBA191126_9J_2-A1-AD001',
#             'CEMBA191126_9J_1-CEMBA191126_9J_2-A1-AD002',
#             'CEMBA191126_9J_1-CEMBA191126_9J_2-A1-AD004',
#             'CEMBA191126_9J_1-CEMBA191126_9J_2-A1-AD006',
#             'CEMBA191126_9J_1-CEMBA191126_9J_2-A1-AD007']


# Load packages
import cooler
import pathlib

# Read and validate cell ids
# cell_list = cooler.fileops.list_coolers(input_scool)
# scool_cell_ids = set([i.split('/')[-1] for i in cell_list])
# for cell_id in cell_ids:
#     if cell_id not in scool_cell_ids:
#         raise KeyError(f'{cell_id} not in the scool file {input_scool}')

# Get chromosome list
_cell_url = f'{input_scool}::/cells/{cell_ids[0]}'
_cell_cool = cooler.Cooler(_cell_url)
chromosomes = _cell_cool.chromnames

print(len(cell_ids), 'cells to process')
print(len(chromosomes), 'chromosomes in each cell.')

# Final targets
rule summary:
    input:
        expand('{cell_id}.cool', cell_id = cell_ids)
    shell:
        'touch Success && rm -rf impute_*_tmp'

# Impute each chromosome of each cell
rule impute_chrom:
    input:
        input_scool
    output:
        temp('impute_{cell_id}_tmp/{chrom}.hdf')
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


# Aggregate chromosome HDF files for the same cells
rule agg_cell:
    input:
        expand('impute_{{cell_id}}_tmp/{chrom}.hdf', chrom=chromosomes)
    output:
        '{cell_id}.cool'
    threads:
        1
    shell:
        'hic-internal aggregate-chromosomes '
        '--chrom_size_path {chrom_size_path} '
        '--resolution {resolution} '
        '--input_dir impute_{wildcards.cell_id}_tmp '
        '--output_path {output}'
