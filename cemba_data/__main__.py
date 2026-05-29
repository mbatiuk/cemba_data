"""
When adding new function:
1. add a func_register_subparser function to register the subparser
2. add a condition in main func about this new func name, import the real func as func in main
"""

import argparse
import logging
import os
import sys
import cemba_data
from cemba_data import __version__

os.environ["NUMEXPR_MAX_THREADS"] = "8"

log = logging.getLogger()

DESCRIPTION = """
YAP (Yet Another Pipeline) is a mapping pipeline for multiple
snmC-seq based single-cell sequencing technologies.

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


def qsub_register_subparser(subparser):
    parser = subparser.add_parser('qsub',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="General qsub helper, need to prepare a command dict file")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--command_file_path",
        type=str,
        required=True,
        nargs='+',
        help="One or space-separated paths of the command file. Accept 2 different command file format: "
             "1. each line in the file is a full command, will be submitted as a single job. "
             "Qsub parameters can only be specified by qsub_global_parms in this way."
             "2. JSON format, a list of word, where each dict is a command (required key \"command\") "
             "and other optional key specify qsub parameters per job."
    )

    parser_req.add_argument(
        "--working_dir",
        type=str,
        required=True,
        help="Working directory of the qsub project"
    )

    parser_req.add_argument(
        "--project_name",
        type=str,
        required=True,
        help="Name of the qsub project"
    )

    parser.add_argument(
        "--wait_until",
        type=str,
        nargs='+',
        required=False,
        help="If provided with a space separate job id(s), "
             "this job will wait until those job finish first."
    )

    parser.add_argument(
        "--total_cpu",
        type=int,
        required=False,
        default=30,
        help="Total concurrent CPU in qsub list, together with total_mem, "
             "determines how many jobs running in parallel."
    )

    parser.add_argument(
        "--total_mem",
        type=int,
        required=False,
        default=9999,
        help="Total concurrent MEM (GBs) in qsub list, together with total_cpu, "
             "determines how many jobs running in parallel. "
             "Default is 9999, which means only use total_cpu to determine."
    )

    parser.add_argument(
        "--qsub_global_parms",
        type=str,
        required=False,
        default='',
        help="Other global qsub parameters you want to pass to each job's qsub script. "
             "This will cover command.json if set repeatedly."
             "These additional parameters should form in one ';' separated string like this: "
             "'-q=queue_name;-cwd;-wd=/path/to/working/dir;-pe smp=20;-l h_vmem=5G'"
    )
    return


def sbatch_register_subparser(subparser):
    parser = subparser.add_parser('sbatch',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Sbatch wrapper")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--project_name",
        type=str,
        required=True,
        help="Name of the qsub project"
    )

    parser_req.add_argument(
        "--command_file_path",
        type=str,
        required=True,
        help="Each row is a command."
    )

    parser_req.add_argument(
        "--working_dir",
        type=str,
        required=True,
        help="Working directory of the qsub project"
    )

    parser_req.add_argument(
        "--time_str",
        type=str,
        required=True,
        help="Time estimate (upper limit) of Sbatch Jobs."
    )

    parser_req.add_argument(
        "--qos",
        type=str,
        default='serial',
        help="QOS of Slurm. "
    )

    parser_req.add_argument(
        "--mem",
        type=str,
        default='300G',
        help="Memory limit for each job."
    )

    parser_req.add_argument(
        "--cpus",
        type=int,
        default=62,
        help="CPUs per task for each job."
    )

    parser_req.add_argument(
        "--email",
        type=str,
        default=None,
        help="If you want to get notice about job status, put your email here."
    )

    parser_req.add_argument(
        "--email_type",
        type=str,
        default='fail',
        help="Type of status you want to get notice, default is fail, means only get notice when job failed."
    )

    parser_req.add_argument(
        "--template",
        type=str,
        choices=['yap'],
        default='yap',
        help="Type of sbatch template to use."
    )

    parser_req.add_argument(
        "--conda_base",
        type=str,
        default='mamba',
        help="Specify local conda installation type."
             "Accepted values are: 'mamba', 'mambaforge', 'conda', 'miniconda', 'anaconda', 'miniforge', 'miniforge3'"
             "If loaded as a module on HPC specify 'module <module_name>' to use for module load <module_name>"
    )

    parser.add_argument(
        "--max_jobs",
        type=int,
        required=False,
        default=None,
        help="Max number of jobs for the same user (not for the same sbatch command). "
    )

    parser.add_argument(
        "--dry_run",
        required=False,
        action='store_true',
        help="Prepare scripts for sbatch manual submission or debug."
    )
    return


def print_default_config_register_subparser(subparser):
    parser = subparser.add_parser('default-mapping-config',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Print out default config of mapping pipeline")
    from .config import MAPPING_MODE_CHOICES
    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        choices=MAPPING_MODE_CHOICES,
        help="Library mode"
    )

    parser.add_argument(
        "--hisat3n_dna_ref",
        type=str,
        required=True,
        help="Path to the hisat-3n DNA reference"
    )

    parser.add_argument(
        "--hisat3n_rna_ref",
        type=str,
        required=False,
        help="Path to the hisat-3n RNA reference"
    )


    parser.add_argument(
        "--genome_fasta",
        type=str,
        required=True,
        help="Path to the genome fasta file"
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=False,
        help="[mct only] Path to the GTF annotation file, required if mode is mct"
    )

    parser.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help="Path to the chrom size file. "
             "Only chromosomes in this file will be considered."
    )

    parser.add_argument(
        "--nome",
        dest='nome',
        action='store_true',
        help='Does this library have NOMe treatment?'
    )

    parser.add_argument(
        "--annotation_path",
        type=str,
        required=False,
        help="[-mhap only] Path to the *_allc.gz"
    )

    parser.set_defaults(nome=False)
    return


def mapping_cell_fastq_register_subparser(subparser):
    parser = subparser.add_parser('mapping-cell-fastq',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Start mapping pipeline from cell FASTQ files.")

    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument(
        "--output_dir",
        '-o',
        type=str,
        required=True,
        help="Output directory of the mapping pipeline"
    )

    parser_req.add_argument(
        "--config_path",
        "-config",
        type=str,
        required=True,
        help="Path to the mapping config, see 'yap default-mapping-config' about how to generate this file."
    )

    parser_req.add_argument(
        "--fastq_pattern",
        "-fq",
        type=str,
        required=True,
        help="Path pattern with wildcard to match all cell-level FASTQ files, pattern with wildcard must be quoted."
    )

    parser_req.add_argument(
        "--cells_per_group",
        type=int,
        required=False,
        default=64,
        help="Number of cells per mapping group directory (Group0, Group1, ...). "
             "Total number of groups = ceil(n_cells / cells_per_group)."
    )
    parser_req.add_argument(
        "--n_jobs",
        type=int,
        required=False,
        default=64,
        help="Number of jobs to run in parallel in each group"
    )
    parser_req.add_argument(
        "--total_memory_gb",
        type=int,
        required=False,
        default=128,
        help="How many GB memory the computer has in total"
    )
    parser.add_argument(
        "--qos",
        type=str,
        default='serial',
        help="QOS parameter for sbatch script"
    )
    parser.add_argument(
        "--conda_base",
        type=str,
        default='mamba',
        help="Conda base"
    )
    return


def demultiplex_register_subparser(subparser):
    parser = subparser.add_parser('demultiplex',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Run demultiplex.")

    parser.add_argument("--fq_dir", type=str, default="fastq", help="Input fastq directory")
    parser.add_argument("--output_dir", "-o", type=str, default="test", help="Output directory")
    parser.add_argument("--n_jobs", type=int, default=16, help="Number of jobs to run in parallel")
    parser.add_argument("--cells_per_group", type=int, default=64,
                        help="Number of cells per batch subdirectory within each plate. "
                             "After demultiplex, cells are distributed into {plate}Group1, {plate}Group2, ... "
                             "subdirs of this size for parallel mapping on HPC.")
    parser.add_argument("--print_only", action="store_true", help="Print commands only")
    parser.add_argument("--config_path", type=str, default=None,
                        help="Path to mapping config .ini file. "
                             "If provided, total_read_pairs_min and total_read_pairs_max "
                             "are read from it.")
    return


def mapping_register_subparser(subparser):
    parser = subparser.add_parser('mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Run mapping pipeline on local machine.")

    parser.add_argument("--output_dir", "-o", type=str, required=True, help="Output directory")
    parser.add_argument("--config_path", type=str, default="mapping_config.ini", help="Path to the mapping config")
    parser.add_argument("--n_jobs", type=int, default=64, help="Number of jobs per snakemake instance")
    parser.add_argument("--total_jobs", type=int, default=12, help="Total number of jobs to run")
    parser.add_argument("--total_memory_gb", type=int, default=None, help="Total memory in GB")
    parser.add_argument("--print_only", action="store_true", help="Print commands only")
    parser.add_argument("--snakemake_template", type=str, default=None, help="Path to the snakemake template")
    parser.add_argument("--qos", type=str, default='serial', help="QOS parameter for sbatch script")
    parser.add_argument("--conda_base", type=str, default='mamba', help="Conda base")
    parser.add_argument("--time_limit", type=str, default='2-00:00:00',
                        help="Wall time limit for each sbatch job (format: D-HH:MM:SS). Default: 2-00:00:00")
    return



def summary_register_subparser(subparser):
    parser = subparser.add_parser('summary',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Mapping summary after the pipeline finished.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--output_dir",
        '-o',
        type=str,
        required=True,
        help="Pipeline output directory."
    )

    parser_req.add_argument(
        "--notebook",
        '-nb',
        type=str,
        required=False,
        default=None,
        help="Notebook template for mapping summary, if not provided, will use yap default template."
    )

    parser_req.add_argument(
        "--mode",
        '-m',
        type=str,
        required=False,
        default=None,
        help="Mapping mode (e.g. mc, mct, or m3c). If not provided, auto-detected from output_dir/mapping_config.ini."
    )

    parser_req.add_argument(
        "--kernel_name",
        '-k',
        type=str,
        required=False,
        default='python3',
        help="kernel name for jupyter notebook (default is python3)"
    )
    return



def link_fastq_register_subparser(subparser):
    parser = subparser.add_parser('link-fastq',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Standardize FASTQ file names using symlinks.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "-i", "--in_fq_dir",
        type=str,
        required=True,
        help="Input directory containing raw FASTQ files."
    )

    parser_req.add_argument(
        "-o", "--out_fq_dir",
        type=str,
        required=True,
        help="Output directory for standardized FASTQ symlinks."
    )

    parser_req.add_argument(
        "-p", "--plate_pattern",
        type=str,
        required=True,
        help="Regex with a capturing group () to extract plate name from filename. "
             "Example: '([a-zA-Z0-9]+)'"
    )

    parser.add_argument(
        "-r", "--read_type_pattern",
        type=str,
        default='(R[12])',
        help="Regex to extract read type (R1/R2) from filename. "
             "The captured value must contain '1' or '2'. "
             "Default '(R[12])' works for standard Illumina files."
    )

    parser.add_argument(
        "-g", "--multiplex_group_pattern",
        type=str,
        default=None,
        help="Legacy option. Regex to extract multiplex group (1-6) from filename. "
             "When provided, group is prepended as '{group}-{plate}' in the symlink name "
             "so that files from different groups remain uniquely named. "
             "The group prefix is stripped automatically during demultiplex; "
             "cell-level files are always named {plate}-{index}-R1.fq.gz."
    )

    parser.add_argument(
        "--lane_pattern",
        type=str,
        default=None,
        help=r"Regex to extract lane from filename (e.g. 'L\d+'). "
             "Ignored if --lane is provided. "
             "If neither is provided, L001 is used as fallback."
    )
    parser.add_argument(
        "-l", "--lane",
        type=str,
        default=None,
        help="Manually assign this lane name to ALL files in in_fq_dir. "
             "Takes precedence over --lane_pattern. "
             "Useful when running link-fastq separately per sequencing run "
             "to differentiate runs as L001, L002, etc."
    )

    parser.add_argument(
        "--recursive",
        action='store_true',
        help="Search for FASTQ files recursively in in_fq_dir."
    )
    return


def link_cell_fastq_register_subparser(subparser):
    parser = subparser.add_parser('link-cell-fastq',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Create symlinks for cell-level FASTQ files.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "-i", "--in_fq_dir",
        type=str,
        required=True,
        help="Input directory containing cell-level FASTQ files."
    )

    parser_req.add_argument(
        "-o", "--out_fq_dir",
        type=str,
        required=True,
        help="Output directory for renamed symlinks."
    )

    parser_req.add_argument(
        "-p", "--plate_pattern",
        type=str,
        required=True,
        help="Regex with one capture group matching the plate ID. "
             "Example: '^([^-]+)'"
    )

    parser_req.add_argument(
        "-w", "--well_pattern",
        type=str,
        required=True,
        help="Regex with one capture group matching the well/index ID. "
             "Example: '-([A-P][0-9]+)-R[12]'"
    )

    parser.add_argument(
        "-r", "--read_type_pattern",
        type=str,
        default='(R[12])',
        help="Regex to extract read type (R1/R2) from filename."
    )

    parser.add_argument(
        "--recursive",
        action='store_true',
        help="Search for FASTQ files recursively in in_fq_dir."
    )
    return


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", help="Show version number and exit",
                        version=__version__)
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar=""
    )

    # add subparsers
    qsub_register_subparser(subparsers)
    sbatch_register_subparser(subparsers)
    print_default_config_register_subparser(subparsers)
    mapping_cell_fastq_register_subparser(subparsers)
    demultiplex_register_subparser(subparsers)
    mapping_register_subparser(subparsers)
    summary_register_subparser(subparsers)
    link_fastq_register_subparser(subparsers)
    link_cell_fastq_register_subparser(subparsers)
    # initiate
    args = None
    if len(sys.argv) > 1:
        # print out version
        if sys.argv[1] in ['-v', '--version']:
            print(cemba_data.__version__)
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
    cur_command = args_vars.pop('command')
    # Do real import here:
    if cur_command == 'qsub':
        from cemba_data.mapping.hpc.qsub import qsub as func
    elif cur_command == 'sbatch':
        from cemba_data.mapping.hpc.sbatch import sbatch_submitter as func
    elif cur_command == 'default-mapping-config':
        from .config import print_default_mapping_config as func
    elif cur_command == 'mapping-cell-fastq':
        from .mapping import mapping_cell_fastq as func
    elif cur_command == 'demultiplex':
        from .demultiplex import demultiplex as func
    elif cur_command == 'mapping':
        from .mapping import mapping as func
    elif cur_command == 'summary':
        from cemba_data.summary import final_summary as func
    elif cur_command == 'link-fastq':
        from cemba_data.link_fastq import link_fastq as func
    elif cur_command == 'link-cell-fastq':
        from cemba_data.link_fastq import link_cell_fastq as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    func(**args_vars)
    return


if __name__ == '__main__':
    main()
