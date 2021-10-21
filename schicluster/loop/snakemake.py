import pathlib
import pandas as pd
import schicluster
from .merge_raw_matrix import make_raw_matrix_cell_table

PACKAGE_DIR = pathlib.Path(schicluster.__path__[0])

with open(PACKAGE_DIR / 'loop/generate_matrix_chunk.Snakefile') as tmp:
    GENERATE_MATRIX_CHUNK_TEMPLATE = tmp.read()
with open(PACKAGE_DIR / 'loop/generate_matrix_group.Snakefile') as tmp:
    GENERATE_MATRIX_SCOOL_TEMPLATE = tmp.read()


def prepare_dir(output_dir, chunk_df, dist, cap, pad, gap, resolution,
                min_cutoff, chrom_size_path, keep_cell_matrix, log_e_str):
    output_dir.mkdir(exist_ok=True)
    cell_table_path = str((output_dir / 'cell_table.csv').absolute())
    chunk_df[['cell_url']].to_csv(cell_table_path)
    parameters = dict(dist=dist,
                      cap=cap,
                      pad=pad,
                      gap=gap,
                      resolution=resolution,
                      min_cutoff=min_cutoff,
                      cell_table_path=f'"{cell_table_path}"',
                      chrom_size_path=f'"{chrom_size_path}"',
                      keep_cell_matrix=keep_cell_matrix,
                      log_e_str=f'"{log_e_str}"')
    parameters_str = '\n'.join(f'{k} = {v}'
                               for k, v in parameters.items())

    with open(output_dir / 'Snakefile', 'w') as f:
        f.write(parameters_str + GENERATE_MATRIX_CHUNK_TEMPLATE)
    return


def prepare_loop_snakemake(cell_table_path, output_dir, chrom_size_path, chunk_size=100, dist=10050000,
                           cap=5, pad=5, gap=2, resolution=10000, min_cutoff=1e-6,
                           keep_cell_matrix=False, cpu_per_job=10, log_e=False, raw_resolution_str='10K'):
    cell_table = pd.read_csv(cell_table_path, index_col=0, sep='\t',
                             names=['cell_id', 'cell_url', 'cell_group'])
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(exist_ok=True)

    # a single dir for raw matrix
    cell_table_raw = make_raw_matrix_cell_table(cell_table_path, raw_resolution_str)
    raw_dir = output_dir / 'raw'
    raw_dir.mkdir(exist_ok=True)
    raw_table_path = raw_dir / 'cell_table.tsv'
    cell_table_raw.to_csv(raw_table_path, sep='\t', header=None)
    if raw_resolution_str == '10K':
        raw_resolution = 10000
    else:
        raise NotImplementedError
    merge_raw_cmd = f'hic-internal merge-raw-scool ' \
                    f'--chrom_size_path {chrom_size_path} ' \
                    f'--resolution {raw_resolution} ' \
                    f'--cell_table_path {raw_table_path} ' \
                    f'--output_dir {raw_dir} ' \
                    f'--cpu {cpu_per_job}'

    if log_e:
        log_e_str = '--log_e'
    else:
        log_e_str = ''

    chunk_parameters = dict(
        dist=dist,
        cap=cap,
        pad=pad,
        gap=gap,
        resolution=resolution,
        min_cutoff=min_cutoff,
        chrom_size_path=chrom_size_path,
        keep_cell_matrix=keep_cell_matrix,
        log_e_str=log_e_str
    )

    total_chunk_dirs = []
    group_chunks = {}
    for group, group_df in cell_table.groupby('cell_group'):
        group_chunks[group] = []
        if group_df.shape[0] <= chunk_size:
            this_dir = output_dir / f'{group}_chunk0'
            prepare_dir(this_dir, group_df, **chunk_parameters)
            total_chunk_dirs.append(this_dir)
            group_chunks[group].append(this_dir)
        else:
            group_df['chunk'] = [i // chunk_size for i in range(group_df.shape[0])]
            for chunk, chunk_df in group_df.groupby('chunk'):
                this_dir = output_dir / f'{group}_chunk{chunk}'
                prepare_dir(this_dir, chunk_df, **chunk_parameters)
                total_chunk_dirs.append(this_dir)
                group_chunks[group].append(this_dir)

    with open(output_dir / 'snakemake_cmd_step1.txt', 'w') as f:
        f.write(merge_raw_cmd + '\n')
        for chunk_dir in total_chunk_dirs:
            cmd = f'snakemake -d {chunk_dir} --snakefile {chunk_dir}/Snakefile -j {cpu_per_job}'
            f.write(cmd + '\n')

    # prepare the second step that merge cell chunks into groups
    scool_parameters = dict(
        output_dir=f'"{output_dir}"',
        chrom_size_path=f'"{chrom_size_path}"',
        resolution=resolution
    )
    parameters_str = '\n'.join(f'{k} = {v}'
                               for k, v in scool_parameters.items())
    with open(output_dir / 'Snakefile', 'w') as f:
        f.write(parameters_str + GENERATE_MATRIX_SCOOL_TEMPLATE)
    with open(output_dir / 'snakemake_cmd_step2.txt', 'w') as f:
        cmd = f'snakemake -d {output_dir} --snakefile {output_dir}/Snakefile -j {cpu_per_job}'
        f.write(cmd + '\n')
    return


def check_chunk_dir_finish(output_dir):
    output_dir = pathlib.Path(output_dir)
    flag_path = output_dir / 'chunk_finished'
    if flag_path.exists():
        return

    path_not_finished = []
    for chunk_dir in output_dir.glob('*_chunk*'):
        if not (chunk_dir / 'finish').exists():
            path_not_finished.append(str(chunk_dir))
    if len(path_not_finished) > 0:
        path_str = '\n'.join(path_not_finished)
        raise ValueError('These chunk dirs have not finish successfully:\n' + path_str)
    else:
        with open(flag_path, 'w') as _:
            pass
    return
