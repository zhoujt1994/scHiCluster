import pathlib
import schicluster
import cooler

PACKAGE_DIR = pathlib.Path(schicluster.__path__[0])


def prepare_impute(input_scool,
                   output_dir,
                   chrom_size_path,
                   output_dist,
                   window_size,
                   step_size,
                   resolution,
                   batch_size=100,
                   logscale=False,
                   pad=1,
                   std=1,
                   rp=0.5,
                   tol=0.01,
                   min_cutoff=1e-5,
                   cpu_per_job=10):
    """
    prepare snakemake files and directory structure for cell contacts imputation
    """
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(exist_ok=True)
    input_scool = str(pathlib.Path(input_scool).absolute())

    with open(PACKAGE_DIR / 'impute/impute.Snakefile') as f:
        snake_template = f.read()
    cell_list = cooler.fileops.list_coolers(input_scool)
    scool_cell_ids = [i.split('/')[-1] for i in cell_list]
    if logscale:
        logscale_str = '--logscale'
    else:
        logscale_str = ''

    chunk_dirs = []
    for i, chunk_start in enumerate(range(0, len(scool_cell_ids), batch_size)):
        chunk_dir = output_dir / f'chunk{i}'
        chunk_dir.mkdir(exist_ok=True)

        this_cell_ids = scool_cell_ids[chunk_start:chunk_start + batch_size]
        parameters = dict(
            input_scool=f"'{pathlib.Path(input_scool).absolute()}'",
            chrom_size_path=f"'{pathlib.Path(chrom_size_path).absolute()}'",
            logscale_str=f'"{logscale_str}"',
            pad=pad,
            std=std,
            window_size=int(window_size),
            step_size=int(step_size),
            resolution=int(resolution),
            output_dist=int(output_dist),
            rp=rp,
            tol=tol,
            min_cutoff=min_cutoff,
            cell_ids=str(this_cell_ids)
        )
        parameters_str = '\n'.join([f'{k} = {v}' for k, v in parameters.items()])
        this_snakefile = parameters_str + snake_template
        with open(output_dir / f'chunk{i}/Snakefile', 'w') as f:
            f.write(this_snakefile)
        chunk_dirs.append(chunk_dir)

    with open(output_dir / 'snakemake_cmd.txt', 'w') as f:
        for chunk_dir in chunk_dirs:
            cmd = f'snakemake -d {chunk_dir} --snakefile {chunk_dir}/Snakefile -j {cpu_per_job}'
            f.write(cmd + '\n')
    return
