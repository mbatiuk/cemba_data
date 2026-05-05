import glob
import pathlib
import subprocess
import os
import pandas as pd
from snakemake.io import glob_wildcards
import cemba_data
from .m3c import m3c_config_str
from .mc import mc_config_str
from .mct import mct_config_str
from ...utilities import get_configuration
from ..hpc import prepare_run

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])

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
	except:
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
		snakefile_path = os.path.join(PACKAGE_DIR, f'files/smk/hisat3n/{mode.lower()}.smk')

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

	total_snakefile = cell_id_str + snake_template
	if not os.path.exists(os.path.join(output_dir, 'snakemake')):
		os.makedirs(os.path.join(output_dir, 'snakemake'), exist_ok=True)
	subprocess.run(['touch', os.path.join(output_dir, 'snakemake/hisat3n')], check=True)

	with open(os.path.join(sub_folder, 'Snakefile'), 'w') as f:
		f.write(total_snakefile)
	return


def mapping_cell_fastq(output_dir, fastq_pattern, config_path, n_group=64,
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

	# make symlink of fastq files, using dir structure of demultiplex
	# the cells are randomly grouped though, max group is 64
	groups = min(n_group, fastq_df.shape[0])
	for i, (cell_id, (r1_path, r2_path)) in enumerate(fastq_df.sample(fastq_df.shape[0]).iterrows()):
		group_id = i % groups
		fastq_dir = output_dir / f'Group{group_id}/fastq'
		fastq_dir.mkdir(exist_ok=True, parents=True)

		# make symlinks
		new_r1_path = fastq_dir / r1_path.name
		new_r1_path.symlink_to(r1_path)
		if pd.notna(r2_path):
			new_r2_path = fastq_dir / r2_path.name
			new_r2_path.symlink_to(r2_path)

	for group_id in range(groups):
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
                snakemake_template=None, qos='serial', conda_base='mamba'):
	"""
	Run mapping pipeline on local machine.
	"""
	output_folder = os.path.expanduser(output_dir)
	if not os.path.exists(output_folder):
		os.makedirs(output_folder, exist_ok=True)
	os.system(f"cp {config_path} {output_folder}/mapping_config.ini")
	
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
		            qos=qos, conda_base=conda_base)
	else:
		for cmd in cmds:
			print(f"{cmd}")
			with open(log_path, 'a') as f:
				f.write(cmd + '\n')
			os.system(cmd)
