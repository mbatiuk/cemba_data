import pathlib

import cemba_data
from ..utilities import MAPPING_MODE_CHOICES

# Load defaults
PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def print_default_mapping_config(mode,
								 barcode_version,
								 genome_fasta,
								 hisat3n_dna_ref=None,
								 hisat3n_rna_ref=None,
								 gtf=None,
								 nome=False,
								 chrom_size_path=None,
								 **kwargs):
	mode = mode.lower()
	if mode.split('-')[0] not in MAPPING_MODE_CHOICES:
		raise ValueError(f'Unknown mode {mode}')

	barcode_version = barcode_version.upper()
	if barcode_version not in ['V1', 'V2']:
		raise ValueError(f'Unknown mode {barcode_version}')

	if not hisat3n_rna_ref is None:
		hisat3n_rna_ref = hisat3n_rna_ref
	if not hisat3n_dna_ref is None:
		hisat3n_dna_ref = hisat3n_dna_ref

	if mode.split('-')[0] == 'mct':
		if hisat3n_rna_ref is None:
			raise ValueError('hisat3n_rna_ref is required if mode is mct.')
		if gtf is None:
			raise ValueError('gtf must be provided when mode is mct.')
		gtf = gtf

	if chrom_size_path is None:
		raise ValueError('chrom_size_path must be provided.')
	chrom_size_path = chrom_size_path

	if mode.split('-')[0] == 'm3c':
		pass

	genome_fasta = genome_fasta

	if mode.split('-')[0] == 'mc':
		if nome:
			config_path = PACKAGE_DIR / 'files/default_config/mapping_config_nome.ini'
		else:
			config_path = PACKAGE_DIR / 'files/default_config/mapping_config_mc.ini'
		with open(config_path) as f:
			config_content = f.read()
	elif mode.split('-')[0] == 'mct':
		if nome:
			config_path = PACKAGE_DIR / 'files/default_config/mapping_config_mct-nome.ini'
		else:
			config_path = PACKAGE_DIR / 'files/default_config/mapping_config_mct.ini'
		with open(config_path) as f:
			config_content = f.read()
		if not hisat3n_rna_ref is None:
			config_content = config_content.replace('CHANGE_THIS_TO_YOUR_HISAT3N_RNA_REFERENCE',
													str(hisat3n_rna_ref))
		if not gtf is None:
			config_content = config_content.replace('CHANGE_THIS_TO_YOUR_GENE_ANNOTATION_GTF', str(gtf))
	elif mode.split('-')[0] == 'm3c':
		config_path = PACKAGE_DIR / 'files/default_config/mapping_config_m3c.ini'
		with open(config_path) as f:
			config_content = f.read()
	else:
		raise ValueError(f'Unknown mode {mode}')

	config_content = config_content.replace('CHANGE_THIS_TO_YOUR_MODE', mode)
	config_content = config_content.replace('CHANGE_THIS_TO_YOUR_CHROM_SIZE_PATH', str(chrom_size_path))
	config_content = config_content.replace('USE_CORRECT_BARCODE_VERSION_HERE', barcode_version)
	if not hisat3n_dna_ref is None:
		config_content = config_content.replace('CHANGE_THIS_TO_YOUR_HISAT3N_DNA_REFERENCE', str(hisat3n_dna_ref))
	config_content = config_content.replace('CHANGE_THIS_TO_YOUR_REFERENCE_FASTA', str(genome_fasta))
	print(config_content)
	for key in kwargs:
		if kwargs[key] is None:
			continue
		print(f"{key} = {kwargs[key]}")
	return
