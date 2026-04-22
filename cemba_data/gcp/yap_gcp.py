import os
import pandas as pd
import cemba_data
import glob
import pathlib
from cemba_data.mapping.pipelines import make_all_snakefile, prepare_run
from snakemake.io import glob_wildcards

PACKAGE_DIR = cemba_data.__path__[0]
from cemba_data.demultiplex.fastq_dataframe import _parse_fastq_path
from cemba_data.demultiplex import _parse_index_fasta


def make_fastq_df(fq_dir):
	"""
	Create a DataFrame for FASTQ files in a local directory.
	"""
	print(f"Scanning directory: {fq_dir}")
	input_files = glob.glob(os.path.join(os.path.expanduser(fq_dir), "*.fastq.gz")) + \
	              glob.glob(os.path.join(os.path.expanduser(fq_dir), "*.fq.gz"))

	df = pd.DataFrame([_parse_fastq_path(file) for file in input_files])
	df.fastq_path = df.fastq_path.apply(lambda x: str(x))
	return df


def get_fastq_info(fq_dir, local_outdir="./"):
	if not os.path.exists(local_outdir):
		os.makedirs(local_outdir, exist_ok=True)
	outfile = os.path.join(local_outdir, "fastq_info.txt")
	if os.path.exists(outfile):
		df = pd.read_csv(outfile, sep='\t')
		return df

	df = make_fastq_df(fq_dir)
	df = df.loc[df.read_type.isin(['R1', 'R2'])]
	if df.groupby(['lane', 'read_type'])['uid'].nunique().nunique() != 1:
		print("Warning: number of uids (or fastq files) are different.")
		print(df.groupby(['lane', 'read_type'])['uid'].nunique())

	df = df.loc[df.read_type == 'R1']
	df.rename(columns={'fastq_path': 'R1'}, inplace=True)
	df['R2'] = df.R1.apply(lambda x: x.replace('_R1_', '_R2_'))
	df.sort_values(['uid', 'lane', 'R1'], inplace=True)
	df.to_csv(outfile, sep='\t', index=False)
	return df


def index_name2multiplex_group(x):
	if 'unknow' not in x.lower():
		return ((int(x[1:]) - 1) % 12) // 2 + 1
	else:
		return 'NA'

def get_lanes_info(outdir):
	#  uid={plate}-{multiplex_group}-{primer_name}
	if os.path.exists("lane_info.txt"):
		df1 = pd.read_csv("lane_info.txt", sep='\t')
		df1.fastq_path = df1.fastq_path.apply(lambda x: eval(x))
		return df1

	uids, plates, multiple_groups, primer_names, lanes, index_names, read_types = glob_wildcards(
		os.path.join(outdir, "{uid}/lanes/{plate}-{multiplex_group}-{primer_name}-{lane}-{index_name}-{read_type}.fq.gz"))
	
	if len(uids) == 0:
		print("Run demultiplex.smk first, then run merge_lanes.smk !")
		return None

	df = pd.DataFrame.from_dict({
		'uid': uids,
		'plate': plates,
		'multiplex_group': multiple_groups,
		'primer_name': primer_names,
		'lane': lanes,
		'index_name': index_names,
		'read_type': read_types
	}).drop_duplicates()

	df['fastq_path'] = df.apply(
		lambda row: os.path.join(outdir, row.uid, "lanes", '-'.join(
			row.loc[['uid', 'lane', 'index_name', 'read_type']].map(str).tolist()) + ".fq.gz"), axis=1)

	# V2: if all FASTQs have the same multiplex group in filename (single multiplex group experiment),
	# derive multiplex group from index name.
	# Otherwise multiplex group is already correct in the filename (Ecker lab standard 6-group design).
	if df['multiplex_group'].nunique() == 1:
		df['real_multiplex_group'] = df.index_name.apply(index_name2multiplex_group)
	else:
		df['real_multiplex_group'] = df.multiplex_group.tolist()
	df = df.loc[df.real_multiplex_group != 'NA']

	# new uid (real uid)
	df['uid'] = df.plate.map(str) + '-' + df.real_multiplex_group.map(str) + '-' + df.primer_name.map(str)

	# Put multiple lanes fastq into one list
	df1 = df.loc[:, ['uid', 'index_name', 'read_type',
	                 'fastq_path']].groupby(
		['uid', 'index_name', 'read_type'], as_index=False).agg(lambda x: x.tolist())
	df1.to_csv("lane_info.txt", sep='\t', index=False)
	return df1

def get_random_index(UIDs, local_outdir="./"):
	if not os.path.exists(local_outdir):
		os.makedirs(local_outdir, exist_ok=True)
	outfile = os.path.join(local_outdir, "random_index.txt")
	if os.path.exists(outfile):
		df_index = pd.read_csv(outfile, sep='\t')
		return df_index

	R = []
	for uid in UIDs:
		multiplex_group = uid.split('-')[-2]
		random_index_fa = os.path.join(PACKAGE_DIR, 'files', 'random_index_v2',
		                               f'random_index_v2.multiplex_group_{multiplex_group}.fa')
		index_seq_dict = _parse_index_fasta(random_index_fa)
		index_names = list(index_seq_dict.keys())
		for index_name in index_names:
			for read_type in ['R1', 'R2']:
				R.append([uid, read_type, index_name])

	df_index = pd.DataFrame(R, columns=['uid', 'read_type', 'index_name'])
	df_index.rename(columns={'uid': 'old_uid'}, inplace=True)
	df_index['real_multiplex_group'] = df_index.index_name.apply(index_name2multiplex_group)
	df_index['uid'] = df_index.loc[:, ['old_uid', 'real_multiplex_group']].apply(
		lambda x: '-'.join([x.old_uid.split('-')[0], str(x.real_multiplex_group), x.old_uid.split('-')[-1]]), axis=1)
	df_index['cell_id'] = df_index['uid'] + '-' + df_index['index_name']
	df_index.to_csv(outfile, sep='\t', index=False)
	return df_index


def run_demultiplex(fq_dir="fastq", outdir="test", n_jobs=16, print_only=False):
	"""
	Run demultiplex on local machine.
	"""
	smk1 = os.path.join(PACKAGE_DIR, 'files/smk/demultiplex.smk')
	cmd = f'snakemake -s {smk1} --scheduler greedy --printshellcmds --rerun-incomplete ' \
	      f'--config fq_dir="{fq_dir}" outdir="{outdir}" -j {n_jobs}'

	print(cmd)
	if not print_only:
		os.system(cmd)


def run_mapping(workd,
                config_path="mapping_config.ini",
                n_jobs=64, total_jobs=12, total_memory_gb=None,
                print_only=False,
                snakemake_template=None, qos='serial', conda_base='mamba'):
	"""
	Run mapping pipeline on local machine.
	"""
	output_folder = os.path.expanduser(workd)
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


def start_from_cell_bam(
		indir="bam", bam_pattern="*.hisat3n_dna.all_reads.deduped.bam",
		output_dir="mapping",
		config_path="mapping_config.ini",
		n_jobs=64, print_only=False,
		snakemake_template=None, n_group=64,
		total_memory_gb=128):
	"""
	Start mapping pipeline from BAM files on local filesystem.
	"""
	output_dir = os.path.expanduser(output_dir)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir, exist_ok=True)
	os.system(f"cp {config_path} {output_dir}/mapping_config.ini")
	
	indir = os.path.abspath(os.path.expanduser(indir))
	bam_paths = glob.glob(os.path.join(indir, bam_pattern))
	groups = min(n_group, len(bam_paths))
	for i, bam_path in enumerate(bam_paths):
		group_id = i % groups
		bam_dir = os.path.join(output_dir, f'Group{group_id}/bam')
		os.makedirs(bam_dir, exist_ok=True)
		
		# make symlink of bam files, using dir structure of demultiplex
		new_bam_path = os.path.join(bam_dir, os.path.basename(bam_path))
		if not os.path.exists(new_bam_path):
			os.symlink(bam_path, new_bam_path)
	
	for group_id in range(groups):
		make_all_snakefile(output_dir, f'Group{group_id}',
		                   snakemake_template=snakemake_template,
		                   pattern=f"bam/{bam_pattern}")
	
	print(f"print_only: {print_only}")
	if print_only:
		prepare_run(output_dir, cores_per_job=n_jobs, total_memory_gb=total_memory_gb)
	else:
		# Define the execution parameters
		common_str = f'--default-resources mem_mb=100 --resources mem_mb=50000 --printshellcmds ' \
					 f'--scheduler greedy --rerun-incomplete -j {n_jobs} '

		# Loop through the groups and run them
		for group_id in range(groups):
			group_path = os.path.join(output_dir, f'Group{group_id}')
			cmd = f"snakemake -s {group_path}/Snakefile {common_str} -d {group_path}"
			print(f"Running: {cmd}")
			os.system(cmd)