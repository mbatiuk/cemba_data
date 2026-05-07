# ==================================================
# Import
# ==================================================
import yaml
from cemba_data.mapping import *

# ==================================================
# Preparation
# ==================================================
# read mapping config and put all variables into the locals()
DEFAULT_CONFIG = {
    'hisat3n_repeat_index_type': '',
    'r1_adapter': 'AGATCGGAAGAGCACACGTCTGAAC',
    'r2_adapter': 'AGATCGGAAGAGCGTCGTGTAGGGA',
    'r1_right_cut': 10,
    'r2_right_cut': 10,
    'r1_left_cut': 10,
    'r2_left_cut': 10,
    'min_read_length': 30,
    'num_upstr_bases': 0,
    'num_downstr_bases': 2,
    'compress_level': 5,
    'hisat_threads': 11,
    'hisat3n_threads': 11,
    # the post_mapping_script can be used to generate dataset, run other process etc.
    # it gets executed before the final summary function.
    # the default command is just a placeholder that has no effect
    'post_mapping_script': 'true',
    'feature_type': 'gene',
    'id_type': 'gene_id',
}
REQUIRED_CONFIG = ['hisat_dna_reference', 'hisat_rna_reference', 'gtf_path', 'reference_fasta', 'chrom_size_path']

bam_dir = "bam"
allc_dir = "allc"
hic_dir = "hic"

local_config = read_mapping_config()
DEFAULT_CONFIG.update(local_config)
if 'hisat3n_dna_reference' in DEFAULT_CONFIG:
    DEFAULT_CONFIG['hisat_dna_reference'] = DEFAULT_CONFIG['hisat3n_dna_reference']
if 'hisat3n_rna_reference' in DEFAULT_CONFIG:
    DEFAULT_CONFIG['hisat_rna_reference'] = DEFAULT_CONFIG['hisat3n_rna_reference']

for k, v in DEFAULT_CONFIG.items():
    if k not in config:
        config[k] = v

missing_key = []
for k in REQUIRED_CONFIG:
    if k not in config:
        missing_key.append(k)
if len(missing_key) > 0:
    raise ValueError('Missing required config: {}'.format(missing_key))



mcg_context = 'CGN' if int(config['num_upstr_bases']) == 0 else 'HCGN'
#repeat_index_flag = "--repeat" if config['hisat3n_repeat_index_type'] == 'repeat' else "--no-repeat-index"
repeat_index_flag="--no-repeat-index" #repeat would cause some randomness, get different output (mapping summary) even using the same input and parameters
allc_mcg_dir=f"allc-{mcg_context}"
# print(f"bam_dir: {bam_dir}\n allc_dir: {allc_dir}\n hic_dir: {hic_dir} \n allc_mcg_dir: {allc_mcg_dir}")

for dir in [bam_dir, allc_dir, allc_mcg_dir]:
    if not os.path.exists(dir):
        os.mkdir(dir)

def get_fastq_path():
    return "fastq/{cell_id}-{read_type}.fq.gz"