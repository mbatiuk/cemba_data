"""
When adding new function:
1. add a func_register_subparser function to register the subparser
2. add a condition in main func about this new func name, import the real func as func in main
"""

import argparse
import inspect
import logging
import sys

import cemba_data
from .__main__ import setup_logging

log = logging.getLogger()

DESCRIPTION = """
yap-internal is used for automation, not intend to be used by end user. 
Use yap instead. 
"""

EPILOG = ''


def dss_two_internal_subparser(subparser):
    parser = subparser.add_parser('dss-two',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Run DSS two-group DMR")
    parser_req = parser.add_argument_group("Required inputs")
    parser_req.add_argument("--allc_table_path", type=str, required=True,
                            help="Three columns separated by tab: 1) allc path; 2) sample id 3) group id")
    parser_req.add_argument("--output_dir", type=str, required=True, help="Output directory")
    parser_req.add_argument("--study_name", type=str, required=True, help="Name of the study")
    parser_req.add_argument("--chrom_sizes_path", type=str, required=True,
                            help="Path to genome chrom size file")
    parser.add_argument("--chroms", type=str, nargs='+', default=None,
                        help="Chromosomes to consider.")
    parser.add_argument("--p_threshold", type=float, default=0.001,
                        help="FDR threshold to select sig DML")
    parser.add_argument("--min_cg", type=int, required=True, default=1,
                        help="Minimum CpGs for a DMR")
    parser.add_argument("--min_len", type=int, required=True, default=1,
                        help="Minimum length for a DMR")
    parser.add_argument("--sig_ratio", type=float, required=True, default=0.5,
                        help="Minimum ratio of CpGs that are significant in a DMR.")
    parser.add_argument("--delta", type=float, required=True, default=0.1,
                        help="Methylation delta that considered to be informative.")
    parser.add_argument("--cpu", type=int, required=True, default=10,
                        help="Number of CPUs to use")
    parser.add_argument("--chunk_size", type=int, default=50000000,
                        help="chunk size to parallel jobs")
    parser.add_argument("--not_smoothing", dest='smoothing', action='store_false',
                        help='Do not perform smoothing (will perform smoothing by default).')
    parser.set_defaults(smoothing=True)
    parser.add_argument("--save_dml", dest='save_dml', action='store_true',
                        help='Save all the DML files (will delete DML table by default).')
    parser.set_defaults(save_dml=False)
    return


def dss_multi_internal_subparser(subparser):
    parser = subparser.add_parser('dss-multi',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Run DSS multi-group DMR")
    parser_req = parser.add_argument_group("Required inputs")
    parser_req.add_argument("--allc_table_path", type=str, required=True,
                            help="Three columns separated by tab: 1) allc path; 2) sample id 3) group id")
    parser_req.add_argument("--output_dir", type=str, required=True, help="Output directory")
    parser_req.add_argument("--study_name", type=str, required=True, help="Name of the study")
    parser_req.add_argument("--chrom_sizes_path", type=str, required=True,
                            help="Path to genome chrom size file")
    parser.add_argument("--chroms", type=str, nargs='+', default=None,
                        help="Chromosomes to consider.")
    parser.add_argument("--p_threshold", type=float, default=0.001,
                        help="FDR threshold to select sig DML")
    parser.add_argument("--min_cg", type=int, required=True, default=1,
                        help="Minimum CpGs for a DMR")
    parser.add_argument("--min_len", type=int, required=True, default=1,
                        help="Minimum length for a DMR")
    parser.add_argument("--sig_ratio", type=float, required=True, default=0.5,
                        help="Minimum ratio of CpGs that are significant in a DMR.")
    parser.add_argument("--cpu", type=int, required=True, default=10,
                        help="Number of CPUs to use")
    parser.add_argument("--chunk_size", type=int, default=50000000,
                        help="chunk size to parallel jobs")
    parser.add_argument("--no_smooth", dest='smoothing', action='store_false',
                        help='Do not perform smoothing (will perform smoothing by default).')
    parser.set_defaults(smoothing=True)
    parser.add_argument("--save_dml", dest='save_dml', action='store_true',
                        help='Save all the DML files (will delete DML table by default).')
    parser.set_defaults(save_dml=False)
    return


def dmrseq_internal_subparser(subparser):
    parser = subparser.add_parser('dmrseq',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Run DMRseq")
    parser_req = parser.add_argument_group("Required inputs")
    parser_req.add_argument("--allc_table_path", type=str, required=True,
                            help="Three columns separated by tab: 1) allc path; 2) sample id 3) group id")
    parser_req.add_argument("--output_dir", type=str, required=True, help="Output directory")
    parser_req.add_argument("--study_name", type=str, required=True, help="Name of the study")
    parser_req.add_argument("--chrom_sizes_path", type=str, required=True,
                            help="Path to genome chrom size file")
    parser.add_argument("--chroms", type=str, nargs='+', default=None,
                        help="Chromosomes to consider.")
    parser.add_argument("--test_covariate", type=str, default='group',
                        help="Test covariate name in the allc path table")
    parser.add_argument("--match_covariate", type=str, default=None, nargs='+',
                        help="Matched covariate name in the allc path table.")
    parser.add_argument("--adjust_covariate", type=str, default=None, nargs='+',
                        help="Adjust covariate name in the allc path table.")
    parser.add_argument("--cutoff", type=float, default=0.1,
                        help="Methylation delta that considered to be informative.")
    parser.add_argument("--min_num_region", type=int, default=3,
                        help="Minimum number of CpGs in a DMR, min is 3.")
    parser.add_argument("--bp_span", type=int, default=1000)
    parser.add_argument("--min_in_span", type=int, default=30)
    parser.add_argument("--max_gap_smooth", type=int, default=2500)
    parser.add_argument("--max_gap", type=int, default=1000)
    parser.add_argument("--max_perms", type=int, default=10)
    parser.add_argument("--stat", type=str, default='stat')
    parser.add_argument("--chrs_per_chunk", type=int, default=1)
    parser.add_argument("--block_size", type=int, default=5000)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--template_path", type=str, default='default')
    parser.add_argument("--chunk_size", type=int, default=50000000,
                        help="chunk size to parallel jobs")
    parser.add_argument("--no_smooth", dest='smooth', action='store_false',
                        help='Do not perform smoothing (will perform smoothing by default).')
    parser.set_defaults(smooth=True)
    parser.add_argument("--block", dest='block', action='store_true',
                        help='Save all the DML files (will delete DML table by default).')
    parser.set_defaults(block=False)
    return


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
    for name, register_subparser_func in inspect.getmembers(current_module, inspect.isfunction):
        if 'internal_subparser' in name:
            register_subparser_func(subparsers)

    # initiate
    if len(sys.argv) > 1:
        if sys.argv[1] in ['-v', '--version']:
            print(cemba_data.__version__)
            exit()
        else:
            args = parser.parse_args()
    else:
        parser.parse_args(["-h"])
        exit()

    # set up logging
    if not logging.root.handlers:
        setup_logging(stdout=True, quiet=False)

    # execute command
    args_vars = vars(args)
    cur_command = args_vars.pop('command')
    if cur_command == 'dss-two':
        from .dmr.dss import run_dss_two_group as func
    elif cur_command == 'dss-multi':
        from .dmr.dss import run_dss_multi_group as func
    elif cur_command == 'dmrseq':
        from .dmr.dmrseq import run_dmrseq as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    log.info(f"# Executing {cur_command}...")
    func(**args_vars)
    log.info(f"# {cur_command} finished.")
    return
