import pathlib
import pandas as pd
from cemba_data.utilities import get_configuration


def write_qsub_commands(output_dir, cores_per_job, total_memory_gb=None,
				script_dir=None):
	if total_memory_gb is None:
		total_memory_gb = 2 * cores_per_job
	config_par = ''
	cmds = {}
	snake_files = list(output_dir.glob('*/Snakefile'))
	for snake_file in snake_files:
		uid = snake_file.parent.name
		cmd = f"""snakemake -d {snake_file.parent} --snakefile {snake_file} {config_par} -j {cores_per_job} --rerun-incomplete --scheduler greedy --default-resources mem_mb=100 \
--resources mem_mb={int(1024 * total_memory_gb)} && rm -rf {snake_file.parent}/.snakemake"""
		cmds[uid] = cmd #--resources mem_mb is the limitation.
	script_path = script_dir / 'snakemake_cmd.txt'
	with open(script_path, 'w') as f:
		try:
			uid_order = pd.read_csv(
				output_dir / 'stats/UIDTotalCellInputReadPairs.csv', index_col=0, header=None
			).squeeze().sort_values(ascending=False)
			for uid in uid_order.index:
				if uid in cmds:
					f.write(cmds.pop(uid) + '\n')
			if len(cmds) != 0:
				print(cmds)
				print(uid_order)
				raise ValueError(f'UIDs in Snakefiles not found in UIDTotalCellInputReadPairs.csv: {list(cmds.keys())}')
		except FileNotFoundError:
			# uid_order file do not exist (when starting from cell FASTQs)
			for cmd in cmds.values():
				f.write(cmd + '\n')
	return script_path


def write_sbatch_commands(output_dir, cores_per_job, script_dir, total_mem_mb, qos):
	outdir = str(output_dir.absolute())
	cmds = {}
	snake_files = list(output_dir.glob('*/Snakefile'))
	for snake_file in snake_files:
		uid = snake_file.parent.name
		cmd = f'snakemake ' \
			  f'-d {outdir}/{snake_file.parent.name} ' \
			  f'--snakefile {outdir}/{snake_file.parent.name}/Snakefile ' \
			  f'-j {cores_per_job} ' \
			  f'--default-resources mem_mb=100 --rerun-incomplete ' \
			  f'--resources mem_mb={total_mem_mb} ' \
			  f'&& test -f "{outdir}/{snake_file.parent.name}/MappingSummary.csv.gz" && rm -rf {outdir}/{snake_file.parent.name}/.snakemake'
		cmds[uid] = cmd
	script_path = script_dir / f'snakemake_{qos}_cmd.txt'
	with open(script_path, 'w') as f:
		try:
			uid_order = pd.read_csv(
				output_dir / 'stats/UIDTotalCellInputReadPairs.csv', index_col=0, header=None
			).squeeze().sort_values(ascending=False)
			for uid in uid_order.index:
				if uid in cmds:
					f.write(cmds.pop(uid) + '\n')
			if len(cmds) != 0:
				print(cmds)
				print(uid_order)
				raise ValueError(f'UIDs in Snakefiles not found in UIDTotalCellInputReadPairs.csv: {list(cmds.keys())}')
		except FileNotFoundError:
			# uid_order file do not exist (when starting from cell FASTQs)
			for cmd in cmds.values():
				f.write(cmd + '\n')
	return script_path


def prepare_qsub(name, snakemake_dir, total_jobs, cores_per_job, total_memory_gb):
	memory_gb_per_core = int(total_memory_gb / cores_per_job) if total_memory_gb is not None else 2
	output_dir = snakemake_dir.parent
	qsub_dir = snakemake_dir / 'qsub'
	qsub_dir.mkdir(exist_ok=True)
	script_path = write_qsub_commands(output_dir, cores_per_job, total_memory_gb,
									  script_dir=qsub_dir)
	qsub_str = f"""
#!/bin/bash
#$ -N yap{name}
#$ -V
#$ -l h_rt=99:99:99
#$ -l s_rt=99:99:99
#$ -wd {qsub_dir}
#$ -e {qsub_dir}/qsub.error.log
#$ -o {qsub_dir}/qsub.output.log
#$ -pe smp 1
#$ -l h_vmem=3G

yap qsub \
--command_file_path {script_path} \
--working_dir {qsub_dir} \
--project_name y{name} \
--total_cpu {int(cores_per_job * total_jobs)} \
--qsub_global_parms "-pe smp={cores_per_job};-l h_vmem={memory_gb_per_core}G"
"""
	qsub_total_path = qsub_dir / 'qsub.sh'
	with open(qsub_total_path, 'w') as f:
		f.write(qsub_str)
	print('#' * 60)
	print(f"IF YOU USE QSUB: ")
	print(f"All snakemake commands need to be executed "
		  f"were included in {qsub_total_path}")
	print(f"You just need to qsub this script to "
		  f"map the whole library in {output_dir}")
	print(f"You can also change the per job parameters in {script_path} "
		  f"or change the global parameters in {qsub_total_path}")
	print(f"Read 'yap qsub -h' if you want to have more options about qsub. "
		  f"Alternatively, you can qsub the commands in {script_path} by yourself, "
		  f"as long as they all get successfully executed.")
	print('#' * 60 + '\n')
	return


def prepare_sbatch(name, snakemake_dir, qos, total_memory_gb=None, cores_per_job=None, conda_base='mamba'):
	input_total_mem_mb = total_memory_gb * 1024 if total_memory_gb is not None else None
	output_dir = snakemake_dir.parent
	sbatch_cores_per_job = cores_per_job if cores_per_job is not None else 62
	time_limit = '2-00:00:00'

	total_mem_mb = input_total_mem_mb if input_total_mem_mb is not None else 400 * 1024
	sbatch_mem = f'{total_mem_mb // 1024}G'

	sbatch_dir = snakemake_dir / 'sbatch'
	sbatch_dir.mkdir(exist_ok=True)

	script_path = write_sbatch_commands(output_dir,
										cores_per_job=sbatch_cores_per_job,
										script_dir=sbatch_dir,
										total_mem_mb=total_mem_mb,
										qos=qos)

	sbatch_cmd = f'yap sbatch ' \
				 f'--project_name {name} ' \
				 f'--command_file_path {script_path} ' \
				 f'--working_dir {sbatch_dir} ' \
				 f'--time_str {time_limit} ' \
				 f'--qos {qos} ' \
				 f'--mem {sbatch_mem} ' \
				 f'--cpus {sbatch_cores_per_job} ' \
				 f'--conda_base "{conda_base}"'
	sbatch_total_path = sbatch_dir / f'sbatch-{qos}-qos.sh'
	with open(sbatch_total_path, 'w') as f:
		f.write(sbatch_cmd)
	print('#' * 60)
	print(f'IF YOU USE SBATCH: ')
	print(f"All snakemake commands need to be executed "
		  f"were included in {sbatch_total_path}")
	print(f"You just need to run this script to "
		  f"map the whole library in {output_dir}. "
		  f"Note that this script will not return until all the mapping job finished. "
		  f"So you better run this script with nohup or in a screen.")
	print(f"You can also change "
		  f"the per job parameters in {script_path} "
		  f"or change the global parameters in {sbatch_total_path}")
	print(f"Read 'yap sbatch -h' if you want to have more options about sbatch. "
		  f"Alternatively, you can submit the job scripts in "
		  f"{sbatch_dir} by yourself, "
		  f"as long as they all get successfully executed.")
	print('#' * 40 + '\n')
	return


def prepare_run(output_dir, total_jobs=12, cores_per_job=10, total_memory_gb=None,
				name=None, qos='serial', conda_base='mamba'):
	output_dir = pathlib.Path(output_dir).absolute()
	config = get_configuration(output_dir / 'mapping_config.ini')
	mode = config['mode']
	if mode.split('-')[0] in ['mc', 'm3c'] and cores_per_job < 4:
		raise ValueError(f'cores must >= 4 to run this pipeline.')
	elif mode.split('-')[0] in ['mct'] and cores_per_job < 10:
		raise ValueError(f'cores must >= 10 to run this pipeline.')

	if name is None:
		name = output_dir.name
	snakemake_dir = output_dir / 'snakemake'
	snakemake_dir.mkdir(exist_ok=True)

	prepare_qsub(name=name,
				 snakemake_dir=snakemake_dir,
				 total_jobs=total_jobs,
				 cores_per_job=cores_per_job,
				 total_memory_gb=total_memory_gb)
	prepare_sbatch(name=name, snakemake_dir=snakemake_dir, qos=qos,
				   total_memory_gb=total_memory_gb, cores_per_job=cores_per_job,
				   conda_base=conda_base)

	print(f"Once all commands are executed successfully, use 'yap summary' to generate final mapping summary.")
	return