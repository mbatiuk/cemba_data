"""
Demultiplex pipeline
"""
import glob
import pathlib
import re
import subprocess
import pandas as pd
import os
from cemba_data.utilities import get_configuration

MODULE_DIR = pathlib.Path(__file__).parent.absolute()


def _parse_fastq_path(path):
    """
    UID pattern: {plate}
    FASTQ name pattern (produced by link_fastq):
      {plate}_{lane}_{read_type}.fastq.gz
    or legacy with multiplex group prefix:
      {group}-{plate}_{lane}_{read_type}.fastq.gz
    plate must not contain '-' characters.
    """
    try:
        basename = os.path.basename(path)[:-len('.fastq.gz')]
        plate_field, lane, read_type = basename.rsplit('_', 2)
        plate = plate_field.split('-')[-1]
        if not re.match(r'L\d+$', lane):
            raise ValueError
    except ValueError:
        raise ValueError(f'Found unknown name pattern in path {path}')
    name_dict = dict(plate=plate,
                     lane=lane,
                     read_type=read_type,
                     fastq_path=path,
                     uid=plate)
    return pd.Series(name_dict)


def _parse_index_fasta(fasta_path):
	records = {}
	with open(fasta_path) as f:
		key_line = True
		for line in f:
			if key_line:
				key = line.lstrip('>').rstrip('\n')
				key_line = False
			else:
				value = line.lstrip('^').rstrip('\n')
				records[key] = value
				key_line = True
	return records


def make_fastq_df(fq_dir):
	"""
	Create a DataFrame for FASTQ files in a local directory.
	"""
	print(f"Scanning directory: {fq_dir}")
	input_files = glob.glob(os.path.join(os.path.expanduser(fq_dir), "*.fastq.gz"))

	#Verify if there are any fastq files
	if not input_files:
		raise FileNotFoundError(f"No .fastq.gz files found in {fq_dir}.")

	# Verify every link to a fastq file is real before proceeding
	for f in input_files:
		if not os.path.exists(f):
			raise FileNotFoundError(f"Broken symbolic link to fastq file: {f}")

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

	# Check that both R1 and R2 are present for every fastq read pair
	pivot = df.groupby(['uid', 'lane'])['read_type'].apply(set)
	missing_r1 = pivot[pivot.apply(lambda s: 'R1' not in s)]
	missing_r2 = pivot[pivot.apply(lambda s: 'R2' not in s)]
	errors = []

	if not missing_r1.empty:
		for (uid, lane) in missing_r1.index:
			r2_path = df.loc[
				(df.uid == uid) & (df.lane == lane) & (df.read_type == 'R2'),
				'fastq_path'
			].values[0]
			expected_r1 = r2_path.replace('_R2.fastq.gz', '_R1.fastq.gz')
			errors.append(
				f"  Found R2:       {r2_path}\n"
				f"  Expected R1 at: {expected_r1}"
			)

	if not missing_r2.empty:
		for (uid, lane) in missing_r2.index:
			r1_path = df.loc[
				(df.uid == uid) & (df.lane == lane) & (df.read_type == 'R1'),
				'fastq_path'
			].values[0]
			expected_r2 = r1_path.replace('_R1.fastq.gz', '_R2.fastq.gz')
			errors.append(
				f"  Found R1:       {r1_path}\n"
				f"  Expected R2 at: {expected_r2}"
			)

	if errors:
		raise FileNotFoundError("Incomplete FASTQ pairs detected:\n" + "\n".join(errors))

	df = df.loc[df.read_type == 'R1']
	df.rename(columns={'fastq_path': 'R1'}, inplace=True)
	df['R2'] = df.R1.str.replace('_R1.fastq.gz', '_R2.fastq.gz', regex=False)

	df.sort_values(['uid', 'lane', 'R1'], inplace=True)
	df.to_csv(outfile, sep='\t', index=False)
	return df


def get_random_index(UIDs, cells_per_group, local_outdir="./"):
	if not os.path.exists(local_outdir):
		os.makedirs(local_outdir, exist_ok=True)
	outfile = os.path.join(local_outdir, "random_index.txt")
	if os.path.exists(outfile):
		df_index = pd.read_csv(outfile, sep='\t')
		return df_index

	random_index_fa = MODULE_DIR / 'random_index_v2.fa'
	index_seq_dict = _parse_index_fasta(random_index_fa)
	index_names = list(index_seq_dict.keys())

	R = []
	for uid in UIDs:
		for index_name in index_names:
			for read_type in ['R1', 'R2']:
				R.append([uid, read_type, index_name])

	df_index = pd.DataFrame(R, columns=['uid', 'read_type', 'index_name'])
	df_index['cell_id'] = df_index['uid'] + '-' + df_index['index_name']
	df_index.to_csv(outfile, sep='\t', index=False)
	return df_index


def demultiplex(fq_dir="fastq", output_dir="out", n_jobs=16, cells_per_group=64,
                config_path=None, print_only=False):
	"""
	Run demultiplex on local machine.
	"""
	smk1 = MODULE_DIR / 'demultiplex.smk'
	cmd = f'snakemake -s {smk1} --scheduler greedy --printshellcmds --rerun-incomplete ' \
	      f'-j {n_jobs} --config fq_dir="{fq_dir}" outdir="{output_dir}" cells_per_group={cells_per_group}'

	if config_path is not None:
		user_config = get_configuration(config_path)
		# Only append to the command if the specific keys exist in the .ini
		if 'total_read_pairs_min' in user_config:
			cmd += f' total_read_pairs_min={int(user_config["total_read_pairs_min"])}'
		if 'total_read_pairs_max' in user_config:
			cmd += f' total_read_pairs_max={int(user_config["total_read_pairs_max"])}'

	print(cmd)
	if not print_only:
		subprocess.run(cmd, shell=True, check=True)


def _read_cutadapt_result(stat_path):
	"""
	Parser of cutadapt output
	"""
	with open(stat_path) as f:
		p = re.compile(
			r"Sequence: .+; Type: .+; Length: \d+; Trimmed: \d+ times")
		series = []
		total_pairs = -1
		for line in f:
			if line.startswith('Total read pairs processed'):
				total_pairs = line.split(' ')[-1]
				total_pairs = int(''.join(re.compile(r'\d').findall(total_pairs)))
			m = p.search(line)
			if m is not None:
				result_dict = {}
				for i in m.group().split('; '):
					k, v = i.split(': ')
					result_dict[k] = v
				result_series = pd.Series(result_dict)
				series.append(result_series)
		total_df = pd.DataFrame(series)
		total_df['Trimmed'] = total_df['Trimmed'].apply(
			lambda c: c.split(' ')[0]).astype(int)
		total_df['TotalPair'] = total_pairs
		total_df['Ratio'] = total_df['Trimmed'] / total_pairs
	return total_df

