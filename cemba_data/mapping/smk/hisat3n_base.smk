"""
Snakemake pipeline for hisat-3n mapping of snm3C-seq data

hg38 normal index uses ~9 GB of memory
repeat index will use more memory
"""
import os,sys
import pandas as pd
import yaml
from cemba_data.mapping import *

# ==================================================
# Preparation
# ==================================================

bam_dir = "bam"
allc_dir = "allc"
allc_multi_dir = "allc-multi"
hic_dir = "hic"
mhap_dir = "mhap"

mcg_context = 'CGN' if int(config['num_upstr_bases']) == 0 else 'HCGN'
#repeat_index_flag = "--repeat" if config['hisat3n_repeat_index_type'] == 'repeat' else "--no-repeat-index"
repeat_index_flag="--no-repeat-index" #repeat would cause some randomness, get different output (mapping summary) even using the same input and parameters
allc_mcg_dir=f"allc-{mcg_context}"

for dir in [bam_dir, allc_dir]:
    if not os.path.exists(dir):
        os.mkdir(dir)

def get_fastq_path():
    return "fastq/{cell_id}-{read_type}.fq.gz"