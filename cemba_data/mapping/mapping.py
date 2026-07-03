import glob
import pathlib
import shutil
import subprocess
import os
import pandas as pd
import pysam
from snakemake.io import glob_wildcards
from cemba_data.utilities import get_configuration
from cemba_data.mapping.hpc import prepare_run

MODULE_DIR = pathlib.Path(__file__).parent

def make_all_snakefile(output_dir, subdir=None,
					   snakemake_template=None, pattern="fastq/{cell_id}-R1.fq.gz"):
	"""

	Parameters
	----------
	output_dir :
	subdir :
	snakemake_template :
	pattern : str
		used to get cell_ids
		fastq/{cell_id}-R1.fq.gz, or "CELL_IDS.tsv"

	Returns
	-------

	"""
	try:
		mapping_config_name = [file for file in os.listdir(output_dir) if file.startswith('mapping_config.')][0]
	except (IndexError, FileNotFoundError):
		raise ValueError(f"Could not find mapping_config.* under {output_dir}")
	config = get_configuration(os.path.join(output_dir, mapping_config_name))
	try:
		mode = config['mode']
	except KeyError:
		raise KeyError('mode not found in the config file.')

	if mode.split('-')[0] == 'mc':
		config_str = mc_config_str(config)
	elif mode.split('-')[0] == 'mct':
		config_str = mct_config_str(config)
	elif mode.split('-')[0] == 'm3c':
		config_str = m3c_config_str(config)
	else:
		print(mode)
		raise ValueError(f'Unknown mode {mode}')

	if snakemake_template is not None:
		snakefile_path = os.path.expanduser(snakemake_template)
	else:
		snakefile_path = MODULE_DIR / 'smk' / f'{mode.lower()}.smk'

	with open(snakefile_path) as f:
		snake_template = f.read()

	if subdir is not None:
		sub_folder = os.path.join(output_dir, subdir)
		if not os.path.exists(sub_folder):
			os.makedirs(sub_folder, exist_ok=True)
	else:
		sub_folder = output_dir

	if pattern == 'CELL_IDS':
		cell_ids = pd.read_csv(os.path.join(sub_folder, pattern), sep='\t', index_col=0).index.tolist()
	else:
		cell_ids = glob_wildcards(os.path.join(sub_folder, pattern))[0]

	if len(cell_ids) == 0: # length should be 64
		raise ValueError(f"No cell fastq were identified under {sub_folder}/fastq")
	cell_id_str = f'CELL_IDS = {cell_ids}\n'

	total_snakefile = cell_id_str + config_str + snake_template
	if not os.path.exists(os.path.join(output_dir, 'snakemake')):
		os.makedirs(os.path.join(output_dir, 'snakemake'), exist_ok=True)
	subprocess.run(['touch', os.path.join(output_dir, 'snakemake/hisat3n')], check=True)

	with open(os.path.join(sub_folder, 'Snakefile'), 'w') as f:
		f.write(total_snakefile)
	return


def mapping_cell_fastq(output_dir, fastq_pattern, config_path, cells_per_group=64,
			  n_jobs=64, total_memory_gb=None, qos='serial', conda_base='mamba'):
	output_dir = pathlib.Path(output_dir).absolute()
	if output_dir.exists():
		raise FileExistsError(f'Output dir {output_dir} already exist, please delete it or use another path.')
	output_dir.mkdir()
	subprocess.run(['cp', config_path, f'{output_dir}/mapping_config.ini'], check=True)
	stats_dir = output_dir / 'stats'
	stats_dir.mkdir(exist_ok=True)

	# parse fastq patterns
	fastq_paths = [pathlib.Path(p).absolute() for p in glob.glob(fastq_pattern)]
	r1_records = {}
	r2_records = {}
	for path in fastq_paths:
		*cell_id, suffix = path.name.split('-')
		cell_id = '-'.join(cell_id)
		if suffix == 'R1.fq.gz':
			if cell_id in r1_records:
				raise ValueError(f'Found duplicated cell ID: {cell_id}, '
								 f'File caused this error: {path}')
			r1_records[cell_id] = path
		elif suffix == 'R2.fq.gz':
			if cell_id in r2_records:
				raise ValueError(f'Found duplicated cell ID: {cell_id}, '
								 f'File caused this error: {path}')
			r2_records[cell_id] = path
		else:
			raise ValueError(
				f'Unable to parse read type. Expect file name ends with "-R1.fq.gz" or "-R2.fq.gz", '
				f'File caused this error: {path}'
			)
	fastq_df = pd.DataFrame({'R1Path': r1_records, 'R2Path': r2_records})

	# distribute cells into groups of ~cells_per_group
	n_cells = fastq_df.shape[0]
	n_groups = max(1, (n_cells + cells_per_group - 1) // cells_per_group)
	for i, (cell_id, (r1_path, r2_path)) in enumerate(fastq_df.sample(n_cells).iterrows()):
		group_id = (i % n_groups) + 1  # 1-based to match demultiplex output
		fastq_dir = output_dir / f'Group{group_id}/fastq'
		fastq_dir.mkdir(exist_ok=True, parents=True)

		# make symlinks
		new_r1_path = fastq_dir / r1_path.name
		new_r1_path.symlink_to(r1_path)
		if pd.notna(r2_path):
			new_r2_path = fastq_dir / r2_path.name
			new_r2_path.symlink_to(r2_path)

	for group_id in range(1, n_groups + 1):  # 1-based
		make_all_snakefile(output_dir, subdir=f'Group{group_id}', pattern="fastq/{cell_id}-R1.fq.gz")
	if total_memory_gb is None:
		total_memory_gb = 2 * n_jobs
	prepare_run(output_dir, cores_per_job=n_jobs, total_memory_gb=total_memory_gb,
	            qos=qos, conda_base=conda_base)
	return


def mapping(output_dir,
                config_path="mapping_config.ini",
                n_jobs=64, total_jobs=12, total_memory_gb=None,
                print_only=False,
                snakemake_template=None, qos='serial', conda_base='mamba',
                time_limit='2-00:00:00'):
	"""
	Run mapping pipeline on local machine.
	"""
	output_folder = os.path.expanduser(output_dir)
	if not os.path.exists(output_folder):
		os.makedirs(output_folder, exist_ok=True)
	shutil.copy(config_path, f"{output_folder}/mapping_config.ini")
	
	all_uids, _ = glob_wildcards(os.path.join(output_folder, "{uid}/fastq/{cell_id}-R1.fq.gz"),
	                                        followlinks=True)
	uids = list(set(all_uids))
	common_str = f'--default-resources mem_mb=100 --resources mem_mb=50000 --printshellcmds ' \
	             f'--scheduler greedy --rerun-incomplete -j {n_jobs} '
	log_path = os.path.join(output_folder, "logs.txt")
	pattern = "fastq/{cell_id}-R1.fq.gz"
	
	cmds = []
	for uid in uids:
		make_all_snakefile(output_folder, uid,
		                   snakemake_template=snakemake_template,
		                   pattern=pattern)
		cmd = f"snakemake -s {output_folder}/{uid}/Snakefile {common_str} -d {output_folder}/{uid}"
		cmds.append(cmd)
	
	if print_only:
		if total_memory_gb is None:
			total_memory_gb = 2 * n_jobs
		prepare_run(output_dir=pathlib.Path(output_folder).absolute(),
		            total_jobs=total_jobs,
		            cores_per_job=n_jobs, total_memory_gb=total_memory_gb,
		            qos=qos, conda_base=conda_base, time_limit=time_limit)
	else:
		for cmd in cmds:
			print(f"{cmd}")
			with open(log_path, 'a') as f:
				f.write(cmd + '\n')
			subprocess.run(cmd, shell=True, check=True)


def bam_read_to_fastq_read(read, read_type=None):
	if read_type is None:
		if read.is_read1:
			read_type = '1'
		else:
			read_type = '2'

	fastq_record = f"@{read.qname}_{read_type}\n" \
				   f"{read.query_sequence}\n" \
				   f"+\n" \
				   f"{read.qual}\n"
	return fastq_record


def separate_unique_and_multi_align_reads(in_bam_path,
										  out_unique_path,
										  out_multi_path,
										  out_unmappable_path=None,
										  unmappable_format='auto',
										  mapq_cutoff=10,
										  qlen_cutoff=30,
										  primary_only=True,
										  read_type=None):
	"""
	Separate unique aligned, multi-aligned, and unaligned reads from hisat-3n bam file.

	Parameters
	----------
	in_bam_path
		Path to hisat-3n bam file.
	out_unique_path
		Path to output unique aligned bam file.
	out_multi_path
		Path to output multi-aligned bam file.
	out_unmappable_path
		Path to output unmappable file.
	unmappable_format
		Format of unmappable file, only "bam" and "fastq" supported.
	mapq_cutoff
		MAPQ cutoff for uniquely aligned reads,
		note that for hisat-3n, unique aligned reads always have MAPQ=60
	qlen_cutoff
		read length cutoff for any reads
	primary_only
		If True, only primary alignments (FLAG 256) are considered for multi-aligned reads.
	read_type
		read type, only None, "1" and "2" supported. If the BAM file is paired-end, use None.
	Returns
	-------
	None
	"""
	if out_unmappable_path is not None:
		if unmappable_format == 'auto':
			if out_unmappable_path.endswith('.bam'):
				unmappable_format = 'bam'
			elif out_unmappable_path.endswith('.fastq'):
				unmappable_format = 'fastq'
			else:
				raise ValueError(f'Unmappable format {unmappable_format} not supported.')
		else:
			if unmappable_format not in ['bam', 'fastq']:
				raise ValueError(f'Unmappable format {unmappable_format} not supported.')

	with pysam.AlignmentFile(in_bam_path, index_filename=None) as bam:
		header = bam.header
		with pysam.AlignmentFile(out_unique_path, header=header, mode='wb') as unique_bam, \
				pysam.AlignmentFile(out_multi_path, header=header, mode='wb') as multi_bam:
			if out_unmappable_path is not None:
				if unmappable_format == 'bam':
					unmappable_file = pysam.AlignmentFile(out_unmappable_path, header=header, mode='wb')
				else:
					unmappable_file = open(out_unmappable_path, 'w')
			else:
				unmappable_file = None

			for read in bam:
				# skip reads that are too short
				if read.qlen < qlen_cutoff:
					continue

				if read.mapq > mapq_cutoff:
					unique_bam.write(read)
				elif read.mapq > 0:
					if primary_only and read.is_secondary:
						# skip secondary alignments if primary_only is True,
						# read.is_secondary is True when FLAG contains 256.
						continue
					multi_bam.write(read)
				else:
					# unmappable reads
					if unmappable_file is not None:
						if unmappable_format == 'bam':
							unmappable_file.write(read)
						else:
							unmappable_file.write(bam_read_to_fastq_read(read, read_type=read_type))

			if unmappable_file is not None:
				unmappable_file.close()
	return


def _build_config_str(config, int_parameters, str_parameters, float_parameters=None):
    """Build a typed config dict and serialize it to a Snakemake config string."""
    typed_config = {}

    def _is_none_ish(val):
        return val is None or str(val).strip().lower() in ('none', '')

    for k, default in int_parameters.items():
        if k in config:
            typed_config[k] = int(config[k])
        else:
            if default != 'required':
                typed_config[k] = default
            else:
                raise ValueError(f'Required parameter {k} not found in config. '
                                 f'You can print the newest mapping config template via "yap default-mapping-config".')

    if float_parameters:
        for k, default in float_parameters.items():
            if k in config:
                typed_config[k] = float(config[k])
            else:
                if default != 'required':
                    typed_config[k] = default
                else:
                    raise ValueError(f'Required parameter {k} not found in config. '
                                     f'You can print the newest mapping config template via "yap default-mapping-config".')

    for k, default in str_parameters.items():
        if k in config:
            value = config[k]
            if _is_none_ish(value):
                typed_config[k] = None
            else:
                typed_config[k] = f"'{value}'"
        else:
            if default == 'required':
                raise ValueError(f'Required parameter {k} not found in config. '
                                 f'You can print the newest mapping config template via "yap default-mapping-config".')
            if _is_none_ish(default):
                typed_config[k] = None
            else:
                typed_config[k] = f"'{default}'"

    config_str = "config = {\n"
    for k, v in typed_config.items():
        # When v is the Python object None, f"{v}" produces the literal None (no quotes)
        config_str += f"    '{k}': {v},\n"
    config_str += "}\n"
    return config_str


def m3c_config_str(config):
    """Change the dtype of parameters and make a appropriate string"""
    int_parameters = {
        'overlap': 6,
        'homopolymer_overlap': 30,
        'r1_left_cut': 9,
        'r1_right_cut': 10,
        'r2_left_cut': 10,
        'r2_right_cut': 10,
        'quality_threshold': 20,
        'min_read_length': 30,
        'mapq_threshold': 10,
        'num_upstr_bases': 0,
        'num_downstr_bases': 2,
        'compress_level': 5,
        'min_gap': 2500,
        'hisat3n_threads': 11,
        'optical_duplicate_pixel_distance': 2500,
    }

    str_parameters = {
        'mode': 'm3c',
        'TruSeq2': 'AGATCGGAAGAGCACACGTCTGAAC',
        'TruSeq1_rc': 'AGATCGGAAGAGCGTCGTGTAGGGA',
        'TruSeq1_rc_short': 'AGATCGGAAGAGCGTCGT',
        'TruSeq1': 'TTCCCTACACGACGCTCTTCCGATCT',
        'PolyG': 'G{30}',
        'PolyA': 'A{30}',
        'PolyT': 'T{30}',
        'hisat3n_dna_reference': 'required',
        'hisat3n_repeat_index_type': 'no-repeat',
        'reference_fasta': 'required',
        'chrom_size_path': 'required',
        'annotation_path': None,
		# parameter for mhap pipeline
        'post_mapping_script': 'true',
		# optional script executed after mapping, before final summary.
		# default 'true' is a no-op placeholder. replace with path to your script.
    }

    return _build_config_str(config, int_parameters, str_parameters)


def mc_config_str(config):
    """Change the dtype of parameters and make a appropriate string"""
    int_parameters = {
        'overlap': 6,
        'homopolymer_overlap': 30,
        'r1_left_cut': 9,
        'r1_right_cut': 10,
        'r2_left_cut': 10,
        'r2_right_cut': 10,
        'quality_threshold': 20,
        'min_read_length': 30,
        'mapq_threshold': 10,
        'num_upstr_bases': 0,
        'num_downstr_bases': 2,
        'compress_level': 5,
        'hisat3n_threads': 11,
        'optical_duplicate_pixel_distance': 2500,
    }

    str_parameters = {
        'mode': 'mc',
        'TruSeq2': 'AGATCGGAAGAGCACACGTCTGAAC',
        'TruSeq1_rc': 'AGATCGGAAGAGCGTCGTGTAGGGA',
        'TruSeq1_rc_short': 'AGATCGGAAGAGCGTCGT',
        'TruSeq1': 'TTCCCTACACGACGCTCTTCCGATCT',
        'PolyG': 'G{30}',
        'PolyA': 'A{30}',
        'PolyT': 'T{30}',
        'hisat3n_dna_reference': 'required',
        'hisat3n_repeat_index_type': 'no-repeat',
        'reference_fasta': 'required',
        'chrom_size_path': 'required',
        'annotation_path': None,
		# parameter for mhap pipeline
        'post_mapping_script': 'true',
		# optional script executed after mapping, before final summary.
		# default 'true' is a no-op placeholder. replace with path to your script.
    }

    return _build_config_str(config, int_parameters, str_parameters)


def mct_config_str(config):
    """Change the dtype of parameters and make a appropriate string"""
    int_parameters = {
        'overlap': 6,
        'r1_left_cut': 10,
        'r1_right_cut': 10,
        'r2_left_cut': 10,
        'r2_right_cut': 10,
        'quality_threshold': 20,
        'min_read_length': 30,
        'mapq_threshold': 10,
        'num_upstr_bases': 0,
        'num_downstr_bases': 2,
        'compress_level': 5,
        'dna_cov_min_threshold': 3,
        'rna_cov_min_threshold': 3,
        'hisat3n_threads': 11,
        'optical_duplicate_pixel_distance': 2500,
    }

    float_parameters = {
        'mc_rate_max_threshold': 0.5,
        'mc_rate_min_threshold': 0.9
    }
	
    str_parameters = {
        'mode': 'mct',
        'r1_adapter': 'AGATCGGAAGAGCACACGTCTGAAC',
        'r2_adapter': 'AGATCGGAAGAGCGTCGTGTAGGGA',
        'hisat3n_dna_reference': 'required',
        'hisat3n_rna_reference': 'required',
        'hisat3n_repeat_index_type': 'no-repeat',
        'reference_fasta': 'required',
        'chrom_size_path': 'required',
        'gtf_path': 'required',
        'feature_type': 'gene',
        'id_type': 'gene_id',
        'nome': 'False',
        'annotation_path': None,
		#parameter for mhap pipeline
        'post_mapping_script': 'true',
		# optional script executed after mapping, before final summary.
		# default 'true' is a no-op placeholder. replace with path to your script.
    }

    return _build_config_str(config, int_parameters, str_parameters, float_parameters)

