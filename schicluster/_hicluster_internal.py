import argparse
import inspect
import logging
import sys
from .__main__ import setup_logging

log = logging.getLogger()

DESCRIPTION = """
hic-internal is used for automation, not intend to be used by end user. 
Use hicluster instead. 
"""

EPILOG = ''


def impute_chromosome_internal_subparser(subparser):
    parser = subparser.add_parser('impute-chromosome',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="RWR imputation for one chromosome in one cell")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--scool_url",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--chrom",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--resolution",
        type=int,
        required=True
    )

    parser_req.add_argument(
        "--output_path",
        type=str,
        required=True
    )

    parser.add_argument(
        '--logscale',
        dest='logscale',
        action='store_true')
    parser.set_defaults(logscale=False)

    parser.add_argument(
        "--pad",
        type=int,
        default=1
    )

    parser.add_argument(
        "--std",
        type=int,
        default=1
    )

    parser.add_argument(
        "--rp",
        type=float,
        default=0.5
    )

    parser.add_argument(
        "--tol",
        type=float,
        default=0.01
    )

    parser.add_argument(
        "--window_size",
        type=int,
        default=500000000
    )

    parser.add_argument(
        "--step_size",
        type=int,
        default=10000000
    )

    parser.add_argument(
        "--output_dist",
        type=int,
        default=500000000
    )

    parser.add_argument(
        "--min_cutoff",
        type=float,
        default=0
    )
    return


def aggregate_chromosomes_internal_subparser(subparser):
    parser = subparser.add_parser('aggregate-chromosomes',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Aggregate chromosome HDFs for one cell")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--resolution",
        type=int,
        required=True
    )

    parser_req.add_argument(
        "--input_dir",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--output_path",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--chrom_wildcard",
        type=str,
        default='{chrom}.hdf'
    )


def calculate_loop_matrix_internal_subparser(subparser):
    parser = subparser.add_parser('calculate-loop-matrix',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Calculate Loop Matrix E and T for single chromosome of one cell")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--cell_url",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--chrom",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--resolution",
        type=int,
        required=True
    )

    parser_req.add_argument(
        "--output_prefix",
        type=str,
        required=True
    )

    parser.add_argument(
        "--dist",
        type=int,
        default=10050000
    )

    parser.add_argument(
        "--cap",
        type=int,
        default=5
    )

    parser.add_argument(
        "--pad",
        type=int,
        default=5
    )

    parser.add_argument(
        "--gap",
        type=int,
        default=2
    )

    parser.add_argument(
        "--min_cutoff",
        type=float,
        default=1e-6
    )

    parser.add_argument('--log_e', dest='log_e', action='store_true',
                        help='Normalize E at log scale')
    parser.set_defaults(log_e=False)


def merge_cell_impute_matrix_internal_subparser(subparser):
    parser = subparser.add_parser('merge-cell-impute-matrix',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Merge Loop Matrix E and T from cells to group")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--cell_urls_path",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--chrom",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--output_prefix",
        type=str,
        required=True
    )
    return


def merge_cell_loop_bkg_internal_subparser(subparser):
    parser = subparser.add_parser('merge-loop-matrix',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Merge Loop Matrix E and T from cells to group")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--output_prefix",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--merge_type",
        type=str,
        required=True
    )


def merge_group_chunks_internal_subparser(subparser):
    parser = subparser.add_parser('merge-group-chunks',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Aggregate group chunks cool files to one scool file per group.")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--resolution",
        type=int,
        required=True
    )

    parser_req.add_argument(
        "--group",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True
    )

    parser.add_argument(
        "--matrix_types",
        type=str,
        nargs='+',
        default=['E', 'E2', 'T', 'T2', 'Q'],
        required=False
    )


def merge_raw_scool_internal_subparser(subparser):
    parser = subparser.add_parser('merge-raw-scool',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Merge single cell raw matrix by cluster stored in scool files.")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--resolution",
        type=int,
        required=True
    )

    parser_req.add_argument(
        "--cell_table_path",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=1
    )


def call_loop_internal_subparser(subparser):
    parser = subparser.add_parser('call-loop',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Call loop using group matrix stored in scool file.")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--group_prefix",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--resolution",
        type=int,
        required=True
    )

    parser_req.add_argument(
        "--output_prefix",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--thres_bl",
        type=float,
        default=1.33
    )

    parser_req.add_argument(
        "--thres_donut",
        type=float,
        default=1.33
    )

    parser_req.add_argument(
        "--thres_h",
        type=float,
        default=1.2
    )

    parser_req.add_argument(
        "--thres_v",
        type=float,
        default=1.2
    )

    parser_req.add_argument(
        "--fdr_thres",
        type=float,
        default=0.1
    )

    parser_req.add_argument(
        "--dist_thres",
        type=float,
        default=20000
    )

    parser_req.add_argument(
        "--size_thres",
        type=float,
        default=1
    )


def internal_main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar=""
    )

    # add subparsers
    current_module = sys.modules[__name__]
    # get all functions in parser
    for name, register_subparser_func in inspect.getmembers(current_module, inspect.isfunction):
        if 'internal_subparser' in name:
            register_subparser_func(subparsers)

    # initiate
    args = None
    if len(sys.argv) > 1:
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
        logging.debug(k, v, type(v), sep='\t')

    cur_command = args_vars.pop('command')
    # Do real import here:
    if cur_command == 'impute-chromosome':
        from .impute.impute_chromosome import impute_chromosome as func
    elif cur_command == 'aggregate-chromosomes':
        from .cool.utilities import aggregate_chromosomes as func
    elif cur_command == 'calculate-loop-matrix':
        from .loop.loop_bkg import calculate_chrom_background_normalization as func
    elif cur_command == 'merge-loop-matrix':
        from .loop.merge_cell_to_group import merge_cells_for_single_chromosome as func
    elif cur_command == 'merge-group-chunks':
        from .loop.merge_cell_to_group import merge_group_chunks_to_group_cools as func
    elif cur_command == 'merge-cell-impute-matrix':
        from .impute.merge_cell_to_group import merge_cells_for_single_chromosome as func
    elif cur_command == 'call-loop':
        from .loop.loop_calling import call_loops as func
    elif cur_command == 'merge-raw-scool':
        from .loop.merge_raw_matrix import merge_raw_scool_by_cluster as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    func(**args_vars)
    return
