# Description

This is an independent fork of cemba_data to map single cell DNA methylation and multiome data generated using snmc-type technology from J. Ecker lab. Original code is available at https://github.com/DingWB/cemba_data and https://github.com/lhqing/cemba_data . Codebase was significantly trimmed and refactored, bug fixes were introduced. To increase simplicity legacy code focused on V1 barcoding, bismark mapping and cloud integration was dropped. This is **local-only** implementation for **version 2 barcoding** (384 spearate random indexes for each cell on a 384 plate) based on **hisat3n mapping**.

# Installation
## Create environment and install
```shell
conda install -y -n base -c conda-forge mamba
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

mamba env create -f https://raw.githubusercontent.com/mbatiuk/cemba_data/master/yap.yaml

conda activate yap

```
## Install this fork of cemba_data:
```shell
pip install git+https://github.com/mbatiuk/cemba_data
```

## Reinstall cemba_data after update on github
```shell
pip uninstall -y cemba_data && pip install git+https://github.com/mbatiuk/cemba_data
```

## HISAT-3N Installation
HISAT-3N documentation is available here: https://daehwankimlab.github.io/hisat2/hisat-3n/

HISAT-3N can be built from source:
```shell
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
cd hisat-3n
git checkout hisat3n
make
```
### Temporary add the hisat-3n directory to your PATH
```shell
export PATH=$PWD:$PATH
```

### Permanently add hisat-3n to your PATH, depending on your shell either bash or zsh:
```shell
echo 'export PATH=$PWD:$PATH' >> ~/.bashrc
source ~/.bashrc
```
```shell
echo 'export PATH=$PWD:$PATH' >> ~/.zshrc 
source ~/.zshrc
```

## Build hisat-3n index
### Non-repeat index

```bash
hisat-3n-build --base-change C,T genome.fa genome
```
### Repeat index
```shell
hisat-3n-build --base-change T,C --repeat-index genome.fa genome
```

### Repeat HISAT-3N integrated index with splice site information
```shell
hisat-3n-build --base-change C,T --repeat-index --ss genome.ss --exon genome.exon genome.fa genome
```

### Index for RNA
```shell
hisat2-build -p 16 genome.fa genome
```

# Usage


## Generate config.ini

### m3c - DNA methylation+ 3C
```shell
yap default-mapping-config --mode m3c --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" --hisat3n_dna_ref  "~/Ref/mm10/mm10_ucsc_with_chrL" > m3c_config.ini
```

### mc - DNA methylation
```shell
yap default-mapping-config --mode mc --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" --hisat3n_dna_ref  "~/Ref/mm10/mm10_ucsc_with_chrL" > mc_config.ini
```

### mct - DNA methylation + RNA
```shell
yap default-mapping-config --mode mct --hisat3n_dna_ref "~/Ref/mm10/mm10_ucsc_with_chrL" --hisat3n_rna_ref "~/Ref/mm10/mm10_ucsc_with_chrL" --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" --gtf "~/Ref/mm10/annotations/gencode.vM23.annotation.gtf" > mct_config.ini
```

### mc nome - DNA methylation + NOMe
```shell
yap default-mapping-config --mode mc --nome \
    --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" \
    --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" \
    --hisat3n_dna_ref "~/Ref/mm10/mm10_ucsc_with_chrL" > mc_nome_config.ini
```

### mc-multi - multi-mapping reads retained during mc
```shell
yap default-mapping-config --mode mc-multi \
    --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" \
    --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" \
    --hisat3n_dna_ref "~/Ref/mm10/mm10_ucsc_with_chrL" > mc_multi_config.ini
```

### mct-multi - multi-mapping reads retained during mct
```shell
yap default-mapping-config --mode mct-multi \
    --genome "~/Ref/mm10/mm10.fa" \
    --hisat3n_dna_ref "~/Ref/mm10/dna_index" \
    --hisat3n_rna_ref "~/Ref/mm10/rna_index" \
    --gtf "~/Ref/mm10/annotation.gtf" \
    --chrom_size_path "~/Ref/mm10/mm10.sizes" > mct_multi_config.ini
```
**Attention**, --hisat3n_dna_ref expects prefix of each hisat genome file, not the directory path.
If genome index directory ~/genomes/mus/hisat/ contains files hisat.3n.CT.1.ht2 hisat.3n.CT.2.ht2 etc you need to specify --hisat3n_dna_ref "~/genomes/mus/hisat/hisat" in the command.

Note that NOMe variants are invoked by adding the `--nome` flag to the base mode during configuration generation.

During NOMe variant the GCH contain open chromatin information, HCN contain normal methylation information. An additional one base is added in the context column of ALLC, so we can distinguish GpC sites with HpC sites.


### Supported Modalities

| Modality | Configuration Command |
| :--- | :--- |
| **`m3c`** | `yap default-mapping-config --mode m3c ...` |
| **`m3c-multi`** | `yap default-mapping-config --mode m3c-multi ...` |
| **`mc-nome`** | `yap default-mapping-config --mode mc --nome ...` |
| **`mc-multi`** | `yap default-mapping-config --mode mc-multi ...` |
| **`mct-nome`** | `yap default-mapping-config --mode mct --nome ...` |
| **`mct-multi`** | `yap default-mapping-config --mode mct-multi ...` |

## Demultiplex

### Old yap pipeline
```shell
yap demultiplex --fastq_pattern "Pool_Remind1_m3c/*.fastq.gz" -o your_cell_level_directory -j 16 --config_path m3c_config.ini
```

### New yap-gcp
```shell
 yap-gcp run_demultiplex --fq_dir="Pool_Remind1_m3c" --outdir="your_cell_level_directory" --n_jobs=16 --print_only=True
```

## Mapping


### Generate snakemake, sbatch and qsub scripts to start the mapping
```shell
yap-gcp run_mapping --workd="your_cell_level_directory" --config_path="m3c_config.ini" --n_jobs=62 --total_memory_gb=400 --qos="serial" --conda_base="mamba" --print_only=True
```
--n_jobs - amount of parallel jobs and requested cpu cores if run using sbatch

--total_memory_gb - total RAM memory requested, either per HPC node or locally

--qos - QOS option passed to sbatch script if run on HPC. This is HPC dependent

--conda_base - type of conda installation if using sbatch
  Accepted values are "mamba", "mambaforge", "conda", "miniconda", "anaconda", "miniforge", "miniforge3"
  If loaded as a module on HPC specify "module <module_name>", e.g. "module mamba"
  You can also provide custom path to your conda installation e.g. "/custom/path/to/conda.sh"

--print_only=True - writes snakemake, sbatch and qsub scripts into your_cell_level_directory/snakemake/
--print_only=False - run mapping interactively in the current shell

### Run mapping from generated script files
```shell
bash your_cell_level_directory/snakemake/qsub/snakemake_cmd.txt
```
```shell
bash your_cell_level_directory/snakemake/sbatch/sbatch.sh
```

### Run mapping interactively in the current shell
```shell
yap-gcp run_mapping --workd="your_cell_level_directory" --config_path="m3c_config.ini" --n_jobs=14 --total_memory_gb=32 --print_only=False
```

# Workflow of m3c run
<img src="doc/files/snm3c_dag.svg" title="DAG for snm3c" width="800px">