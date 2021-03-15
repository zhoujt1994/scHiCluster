"""
CLI defined here

When adding new function:
1. add a func_register_subparser function to register the subparser
2. add a condition in main func about this new func name, import the real func as func in main
"""
import numpy as np
import argparse
import inspect
import subprocess
import sys
import logging

log = logging.getLogger()

DESCRIPTION = """
scHiCluster is a toolkit for single-cell HiC data preprocessing, imputation, and clustering analysis.

Current Tool List in scHiCluster:

"""

EPILOG = ''


class NiceFormatter(logging.Formatter):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).
    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = '{}: {}'.format(record.levelname, record.msg)
        return super().format(record)


def validate_environment():
    try:
        # TODO add env validation here
        subprocess.run(['bedtools', '--version'],
                       stderr=subprocess.PIPE,
                       stdout=subprocess.PIPE,
                       encoding='utf8',
                       check=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr)
        raise e
    return True


def setup_logging(stdout=False, quiet=False, debug=False):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Attach handler to the global logger object
    """
    # Due to backwards compatibility, logging output is sent to standard output
    # instead of standard error if the -o option is used.
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setFormatter(NiceFormatter())
    # debug overrides quiet
    if debug:
        level = logging.DEBUG
    elif quiet:
        level = logging.ERROR
    else:
        level = logging.INFO
    stream_handler.setLevel(level)
    log.setLevel(level)
    log.addHandler(stream_handler)


def _str_to_bool(v: str) -> bool:
    if v.lower() in {'true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh'}:
        return True
    else:
        return False


def concat_cells_register_subparser(subparser):
    parser = subparser.add_parser('concat-cells',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="")

    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument(
        "--indir",
        type=str,
        required=True,
        help='Directory of imputed matrices end with /'
    )

    parser_req.add_argument(
        "--chrom",
        type=str,
        required=True,
        help='Chromosome to impute'
    )

    parser_req.add_argument(
        '--mode',
        type=str,
        default=True,
        help='Suffix of imputed matrix file names'
    )

    parser_req.add_argument(
        '--res',
        type=int,
        default=None,
        help='Bin size as integer to generate contact matrix'
    )

    parser.add_argument(
        '--dist',
        type=int,
        default=10000000,
        help='Maximum distance threshold of contacts to use'
    )

    parser.add_argument(
        "--skip_raw",
        dest='save_raw',
        action='store_false',
        help='Whether to skip saving cell-by-feature matrix before SVD'
    )
    parser.set_defaults(save_raw=True)


def generate_matrix_register_subparser(subparser):
    parser = subparser.add_parser('generate-matrix',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="")

    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument(
        '--infile',
        type=str,
        default=None,
        help='Path to the short format contact file'
    )

    parser_req.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='Output directory end with /'
    )

    parser_req.add_argument(
        '--cell',
        type=str,
        default=None,
        help='Specific identifier of a cell'
    )

    parser_req.add_argument(
        '--res',
        type=int,
        default=None,
        help='Bin size as integer to generate contact matrix'
    )

    parser_req.add_argument(
        '--genome',
        type=str,
        default=None,
        help='Genome assembly version'
    )  ## mm10 or hg38

    parser_req.add_argument(
        '--dist',
        type=int,
        default=None,
        help='Minimum distance threshold of contacts to use'
    )


def impute_cell_register_subparser(subparser):
    parser = subparser.add_parser('impute-cell',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="")

    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument(
        '--indir', type=str, required=True, default=None, help='Directory of the contact matrix')
    parser_req.add_argument(
        '--outdir', type=str, required=True, default=None, help='Output directory end with /')
    parser_req.add_argument(
        '--cell', type=str, required=True, default=None, help='Specific identifier of a cell')
    parser_req.add_argument(
        '--chrom', type=str, required=True, default=None, help='Chromosome to impute')
    parser_req.add_argument(
        '--res', type=int, required=True, default=None,
        help='Bin size as integer to generate contact matrix')
    parser_req.add_argument(
        '--genome', type=str, required=True, default=None,
        help='Genome assembly version')  ## mm10 or hg38
    parser_req.add_argument(
        '--mode', type=str, required=True, default=None, help='Suffix of output file name')

    parser.add_argument(
        "--logscale", dest='logscale', action='store_true', help='Whether to log transform raw count or not'
    )
    parser.set_defaults(logscale=False)
    parser.add_argument(
        '--pad', type=int, default=1, help='Gaussian kernal size')
    parser.add_argument(
        '--std', type=float, default=1, help='Gaussian kernal standard deviation')
    parser.add_argument(
        '--rp', type=float, default=0.5, help='Restart probability of RWR')
    parser.add_argument(
        '--tol', type=float, default=0.01, help='Convergence tolerance of RWR')
    parser.add_argument(
        '--rwr_dist', type=int, default=500000000,
        help='Maximum distance threshold of contacts when doing RWR')
    parser.add_argument(
        '--rwr_sparsity', type=float, default=1, help='Minimum sparsity to apply rwr_dist')
    parser.add_argument(
        '--output_dist', type=int, default=500000000,
        help='Maximum distance threshold of contacts when writing output file')
    parser.add_argument(
        '--output_format', type=str, default='hdf', help='Output file format (hdf5 or npz)')


def loop_sc_register_subparser(subparser):
    parser = subparser.add_parser('loop-sc',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="")

    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument('--outdir', type=str, required=True, default=None, help='Directory of imputed matrix')
    parser_req.add_argument('--cell', type=str, required=True, default=None, help='Specific identifier of a cell')
    parser_req.add_argument('--chrom', type=str, required=True, default=None, help='Chromosome imputed')
    parser_req.add_argument('--impute_mode', type=str, required=True, default=None,
                            help='Suffix of imputed matrix file names')
    parser_req.add_argument('--res', type=int, required=True, default=None,
                            help='Bin size as integer to generate contact matrix')
    parser.add_argument('--dist', type=int, default=10000000, help='Maximum distance threshold of contacts to use')
    parser.add_argument('--cap', type=int, default=5, help='Trim Z-scores over the threshold')
    parser.add_argument('--pad', type=int, default=5, help='One direction size of larger square for donut background')
    parser.add_argument('--gap', type=int, default=2, help='One direction size of smaller square for donut background')
    parser.add_argument('--norm_mode', type=str, default='dist_trim', help='Suffix of normalized file names')


def merge_cell_register_subparser(subparser):
    parser = subparser.add_parser('merge-cell',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="")

    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument('--indir', type=str, required=True, default=None,
                            help='Directory of imputed matrices end with /')
    parser_req.add_argument('--cell_list', type=str, required=True, default=None,
                            help='List of cell identifiers to be merged')
    parser_req.add_argument('--group', type=str, required=True, default=None, help='Name of cell group to be merged')
    parser_req.add_argument('--chrom', type=str, required=True, default=None, help='Chromosome to impute')
    parser_req.add_argument('--res', type=int, required=True, default=None,
                            help='Bin size as integer to generate contact matrix')
    parser_req.add_argument('--impute_mode', type=str, required=True, default=None,
                            help='Suffix of imputed matrix file names')

    parser.add_argument('--norm_mode', type=str, default='dist_trim', help='Suffix of normalized file names')
    parser.add_argument('--min_dist', type=int, default=50000, help='Minimum distance threshold of loop')
    parser.add_argument('--max_dist', type=int, default=10000000, help='Maximum distance threshold of loop')
    parser.add_argument('--pad', type=int, default=5, help='One direction size of larger square for donut background')
    parser.add_argument('--gap', type=int, default=2, help='One direction size of smaller square for donut background')
    parser.add_argument('--thres_bl', type=int, default=1.33,
                        help='Fold change threshold against bottom left background')
    parser.add_argument('--thres_d', type=int, default=1.33, help='Fold change threshold against donut background')
    parser.add_argument('--thres_h', type=int, default=1.2, help='Fold change threshold against horizontal background')
    parser.add_argument('--thres_v', type=int, default=1.2, help='Fold change threshold against vertical background')
    return


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     )
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar=""
    )

    # add subparsers
    current_module = sys.modules[__name__]
    # get all functions in parser
    for name, register_subparser_func in inspect.getmembers(current_module, inspect.isfunction):
        if 'register_subparser' in name:
            register_subparser_func(subparsers)

    # initiate
    args = None
    if len(sys.argv) > 1:
        # print out version
        if sys.argv[1] in ['-v', '--version']:
            # print(ALLCools.__version__)
            exit()
        else:
            args = parser.parse_args()
    else:
        # print out help
        parser.parse_args(["-h"])
        exit()

    # set up logging
    if not logging.root.handlers:
        setup_logging(stdout=True,
                      quiet=False)
    # execute command
    args_vars = vars(args)
    for k, v in args_vars.items():
        log.info(f'{k}\t{v}')

    cur_command = args_vars.pop('command').lower().replace('_', '-')
    # Do real import here:
    if cur_command in ['concat-cells']:
        from .dev.concat_cell import concat_cell as func
    elif cur_command in ['generate-matrix']:
        from .dev.generate_matrix import generate_matrix as func
    elif cur_command in ['impute-cell']:
        from .dev.imputecell import impute_cell as func
    elif cur_command in ['loop-sc']:
        from .dev.loop_sc import loop_sc as func
    elif cur_command in ['merge-cell']:
        from .dev.merge_cell import merge_cell as func
    else:
        log.debug(f'{cur_command} is not an valid sub-command')
        parser.parse_args(["-h"])
        return

    validate_environment()

    # run the command
    log.info(f"hicluster: Executing {cur_command}...")
    func(**args_vars)
    log.info(f"{cur_command} finished.")
    return


if __name__ == '__main__':
    main()
