# Description

This is an independent fork of cemba_data to map single cell DNA methylation and multiome data generated using snmc-type technology from J. Ecker lab. Original code is available at https://github.com/DingWB/cemba_data and https://github.com/lhqing/cemba_data . Codebase was significantly trimmed and refactored, bug fixes were introduced. To increase simplicity legacy code focused on V1 barcoding, Index primer, bismark mapping and cloud integration, multiplex groups were dropped. This is **local/HPC-only** implementation for **version 2 barcoding** (384 separate random indexes for each cell on a 384 plate) based on **hisat3n mapping**. In case you need to use multiplex groups - check legacy parameter `--multiplex_group_pattern` in `link-fastq` function. ALLCools dependency was also dropped, while bam-to-allc and extract-allc functions were ported here as a submodule.

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
**Note:** If you use lambda phage DNA as a control for bisulphite conversion you can append lambda DNA fasta at the bottom of genome fasta file. Name lambda DNA as chrL.

### Non-repeat index

```bash
hisat-3n-build --base-change C,T genome.fa hisat
```
### Repeat index
```shell
hisat-3n-build --base-change C,T --repeat-index genome.fa hisat
```
**Note:** repeat index requires more RAM and brings stochasticity into mapping, so results will vary between mapping re-runs.

### Repeat HISAT-3N integrated index with splice site information
```shell
hisat-3n-build --base-change C,T --repeat-index \
    --ss genome.ss --exon genome.exon genome.fa hisat
```

### Index for RNA
```shell
hisat2-build -p 16 genome.fa hisat_rna
```

# Usage


## Generate config.ini

### m3c - DNA methylation+ 3C
```shell
yap default-mapping-config --mode m3c \
    --genome_fasta "~/genomes/mus/mm10_chrL.fa" \
    --chrom_size_path "~/genomes/mus/mm10.sizes" \
    --hisat3n_dna_ref  "~/genomes/mus/hisat/hisat" > m3c_config.ini
```

### mc - DNA methylation
```shell
yap default-mapping-config --mode mc \
    --genome_fasta "~/genomes/mus/mm10_chrL.fa" \
    --chrom_size_path "~/genomes/mus/mm10.sizes" \
    --hisat3n_dna_ref  "~/genomes/mus/hisat/hisat" > mc_config.ini
```

### mct - DNA methylation + RNA
```shell
yap default-mapping-config --mode mct \
    --hisat3n_dna_ref "~/genomes/mus/hisat/hisat" \
    --hisat3n_rna_ref "~/genomes/mus/hisat_rna/hisat_rna" \
    --genome_fasta "~/genomes/mus/mm10_chrL.fa" \
    --chrom_size_path "~/genomes/mus/mm10.sizes" \
    --gtf "~/genomes/mus/gencode.vM23.annotation.gtf" > mct_config.ini
```

### mc nome - DNA methylation + NOMe
```shell
yap default-mapping-config --mode mc --nome \
    --genome_fasta "~/genomes/mus/mm10_chrL.fa" \
    --chrom_size_path "~/genomes/mus/mm10.sizes" \
    --hisat3n_dna_ref "~/genomes/mus/hisat/hisat" > mc_nome_config.ini
```

### mc-multi - multi-mapping reads retained during mc
```shell
yap default-mapping-config --mode mc-multi \
    --genome_fasta "~/genomes/mus/mm10_chrL.fa" \
    --chrom_size_path "~/genomes/mus/mm10.sizes" \
    --hisat3n_dna_ref "~/genomes/mus/hisat/hisat" > mc_multi_config.ini
```

### mc-multi nome - multi-mapping reads retained during mc + NOME
```shell
yap default-mapping-config --mode mc-multi --nome\
    --genome_fasta "~/genomes/mus/mm10_chrL.fa" \
    --chrom_size_path "~/genomes/mus/mm10.sizes" \
    --hisat3n_dna_ref "~/genomes/mus/hisat/hisat" > mc_multi_config.ini
```

### mct-multi - multi-mapping reads retained during mct
```shell
yap default-mapping-config --mode mct-multi \
    --genome_fasta "~/genomes/mus/mm10_chrL.fa" \
    --hisat3n_dna_ref "~/genomes/mus/hisat/hisat" \
    --hisat3n_rna_ref "~/genomes/mus/hisat_rna/hisat_rna" \
    --gtf "~/genomes/mus/annotation.gtf" \
    --chrom_size_path "~/genomes/mus/mm10.sizes" > mct_multi_config.ini
```

### mct-multi nome - multi-mapping reads retained during mct + NOME
```shell
yap default-mapping-config --mode mct-multi --nome\
    --genome_fasta "~/genomes/mus/mm10_chrL.fa" \
    --hisat3n_dna_ref "~/genomes/mus/hisat/hisat" \
    --hisat3n_rna_ref "~/genomes/mus/hisat_rna/hisat_rna" \
    --gtf "~/genomes/mus/annotation.gtf" \
    --chrom_size_path "~/genomes/mus/mm10.sizes" > mct_multi_config.ini
```

**Attention**, `--hisat3n_dna_ref` expects prefix of each hisat genome file, not the directory path.
If genome index directory `~/genomes/mus/hisat/` contains files `hisat.3n.CT.1.ht2` `hisat.3n.CT.2.ht2` etc you need to specify `--hisat3n_dna_ref "~/genomes/mus/hisat/hisat"` in the command.

Note that NOMe variants are invoked by adding the `--nome` flag to the base mode during configuration generation.

During NOMe variant the GCH contain open chromatin information, HCN contain normal methylation information. An additional one base is added in the context column of ALLC, so we can distinguish GpC sites with HpC sites.


### Supported Modalities

| Modality | Configuration Command |
| :--- | :--- |
| **`m3c`** | `yap default-mapping-config --mode m3c ...` |
| **`m3c-multi`** | `yap default-mapping-config --mode m3c-multi ...` |
| **`mc`** | `yap default-mapping-config --mode mc ...` |
| **`mc nome`** | `yap default-mapping-config --mode mc --nome ...` |
| **`mc-multi`** | `yap default-mapping-config --mode mc-multi ...` |
| **`mc-multi nome`** | `yap default-mapping-config --mode mc-multi --nome ...` |
| **`mct`** | `yap default-mapping-config --mode mct ...` |
| **`mct nome`** | `yap default-mapping-config --mode mct --nome ...` |
| **`mct-multi`** | `yap default-mapping-config --mode mct-multi ...` |
| **`mct-multi nome`** | `yap default-mapping-config --mode mct-multi --nome ...` |


## Prepare FASTQ files by linking them with correct naming convention

Raw sequencer outputs often have complex or inconsistent naming conventions. Use `link-fastq` to create symbolic links with standardized names: `{plate}_{lane}_{read_type}.fastq.gz` (or `{group}-{plate}_{lane}_{read_type}.fastq.gz` when `-g` is provided for legacy multiplex-group workflows).

Make sure original raw fastq file names follow these minimal expectations:
- `plate` name must not contain `-`, `_`, or `.` characters. If present they will be stripped.
- `lane` is absent or contains number. `L1`, `1`, `L001`, `l10`, or even `L10000` will be standardized to `L001`, `L010` or `L10000`
- `read_type` must contain `1` or `2` for correct standardization into `R1` and `R2`


```shell
yap link-fastq \
    --in_fq_dir "raw_fastq_dir" \
    --out_fq_dir "fastq" \
    --plate_pattern "([a-zA-Z0-9]+)" \
    --read_type_pattern "(R[12])" \
    --multiplex_group_pattern "(\d+)" \
    --lane_pattern "L\d+" \
    --recursive
```

**Parameters:**
- `-i` / `--in_fq_dir`: (**required**) Input directory containing raw FASTQ files.
- `-o` / `--out_fq_dir`: (**required**) Output directory for standardized FASTQ symlinks.
- `-p` / `--plate_pattern`: (**required**) Regex to extract plate name. Use `()` to capture a sub-portion; without `()` the full match is used.
- `-r` / `--read_type_pattern`: Regex for R1/R2. Default: `(R[12])`. The matched value must contain `1` or `2`.
- `-g` / `--multiplex_group_pattern`: (optional, legacy) Regex to extract multiplex group from filename. Allowed values are 1-6 digits. When provided, the group is prepended as `{group}-{plate}` in the symlink name. This creates unique symlink names in case of multiple .fastq files per single plate. After demultiplexing group is dropped from cell level .fastq file name. 
- `--lane_pattern`: (optional) Regex to extract lane from filename. Use `()` to capture a sub-portion; without `()` the full match is used. Ignored if `--lane` is set.
- `-l` / `--lane`: (optional) Manually assign a lane name to all files in `in_fq_dir`. Number should be present in the `--lane` argument, e.g. `--lane 1` or `--lane L001`. If you have files from multiple sequencing runs of the same library pool you can call `link-fastq` separately per run and use `--lane` to assign a distinct lane number to each run. The pipeline will then merge fastq files correctly. `L001` is used when neither `--lane_pattern` nor `--lane` is provided.
- `--recursive`: Search `in_fq_dir` recursively for FASTQ files.

### How to write patterns (Regular Expressions)

The tool uses Regular Expressions (regex) to find metadata in filenames. If the pattern contains parentheses `()`, only the part inside is extracted. Without `()`, the entire match is used as the value.

For example, to extract the lane from `Sample_L001_R1.fq.gz`:
- `--lane_pattern "_L\d+_"` — matches `_L001_`; the full match including underscores is used — **wrong**.
- `--lane_pattern "_(L\d+)_"` — underscores are context; only `L001` inside `()` is extracted — **correct**.

Use `()` when you need to strip surrounding delimiters or context from the value.

### Common examples of patterns

| Field | Pattern | Extracted value |
| :--- | :--- | :--- |
| **Plate** | `^([a-zA-Z0-9]+)` | Everything alphanumeric at start of filename |
| **Plate** | `([^-_]+)` | Everything until first `-` or `_` |
| **Read Type** | `(R[12])` | `R1` or `R2` (default, usually no need to change) |
| **Lane** | `(L\d+)` | `L001`, `L002`, etc. |
| **Multiplex Group** | `-(\d+)-` | Number between two hyphens |

### Regex quick reference

| Symbol | Meaning | Example |
| :--- | :--- | :--- |
| `.` | Any single character (except newline) | `P.ate` matches `Plate` or `P1ate` |
| `\d` | Any digit (0-9) | `\d+` matches `1` in `123` |
| `\w` | Any word character (letter, digit, or `_`) | `\w` matches `P` in `Plate` |
| `\s` | Any whitespace character | rarely needed in filenames |
| `[abc]` | Any one of the listed characters | `[abc]` matches `a`, `b`, or `c` |
| `[a-zA-Z0-9]` | Any letter or digit | matches `P` in `Plate` or `1` in `123` |
| `[^_-]` | Any character **except** those listed | `[^_-]` captures any character except `-` or `_` |
| `^` | Start of string | `^025` matches `025` only at start |
| `$` | End of string | `\.fastq\.gz$` matches only at end |
| `\|` | Alternation — match either side | `R1\|R2` matches `R1` or `R2` |
| `*` | Zero or more of the preceding | `\d*` matches `` or `123` |
| `+` | One or more of the preceding | `\d+` matches `1` and `123`, not empty |
| `?` | Zero or one of the preceding (optional) | `colou?r` matches `color` or `colour` |
| `{n}` | Exactly n repetitions | `\d{3}` matches exactly `001` |
| `{n,m}` | Between n and m repetitions | `\d{1,3}` matches `1`, `12`, or `123` |
| `()` | Capturing group — **only this part is extracted** | `_(\d+)_` extracts `1` from `_1_` |
| `\.` | Literal dot (escape special characters with `\`) | `\.fastq` matches `.fastq`, not `Xfastq` |

**Note:** if any pattern **starts with `-`**, escape the leading dash with `\`: e.g. `"\-(\d+)-"` instead of `"-(\d+)-"`. This applies to all pattern arguments.


## Prepare cell-level FASTQ files from external sources by linking them with correct naming convention

Use `link-cell-fastq` when working with already-demultiplexed cell-level FASTQ files from an external source where filenames do not follow the canonical format. It creates symlinks renamed to `{plate}-{well}-R1.fq.gz` / `R2.fq.gz` expected by `mapping-cell-fastq`.

For example `025307p2-C2-A13-R1.fq.gz` (external, contains index primer `C2`) becomes `025307p2-A13-R1.fq.gz`.

```shell
yap link-cell-fastq \
    --in_fq_dir "external_cell_fastq" \
    --out_fq_dir "cell_fastq" \
    --plate_pattern "^([^-]+)" \
    --well_pattern "-([A-Z][0-9]+)-R[12]" \
    --read_type_pattern "(R[12])" \
    --recursive
```

**Parameters:**
- `-i` / `--in_fq_dir`: (**required**) Input directory containing cell-level FASTQ files.
- `-o` / `--out_fq_dir`: (**required**) Output directory for renamed symlinks.
- `-p` / `--plate_pattern`: (**required**) Regex with one capture group matching the plate ID.
- `-w` / `--well_pattern`: (**required**) Regex with one capture group matching the well ID.
- `-r` / `--read_type_pattern`: Regex for R1/R2. Default: `(R[12])`. The matched value must contain `1` or `2`.
- `--recursive`: Search `in_fq_dir` recursively for FASTQ files.

The patterns follow the same regex rules described above for `link-fastq`. See the pattern table and regex quick reference above.


## Demultiplex

```shell
yap demultiplex --fq_dir "fastq_dir" \
    --output_dir "your_cell_level_directory" \
    --config_path "m3c_config.ini" \
    --cells_per_group 64 \
    --n_jobs 16 --print_only
```

`--config_path` - path to the `.ini` config file generated by `yap default-mapping-config`. If provided, `total_read_pairs_min` and `total_read_pairs_max` are read from it and passed to the demultiplex pipeline, overriding the built-in defaults (`1` and `10000000`). If omitted, built-in defaults apply.

`--cells_per_group` - number of cells per batch subdirectory. After demultiplex, cells are distributed into `{plate}Group1`, `{plate}Group2`, ... subdirs of this size for parallel mapping on HPC. Default: 64.

`--print_only` - prints demultiplex snakemake command but does not start it interactively

## Mapping


### Generate snakemake, sbatch and qsub scripts to start the mapping
```shell
yap mapping --output_dir "your_cell_level_directory" \
    --config_path "m3c_config.ini" --n_jobs 62 --total_memory_gb 400 \
    --qos "serial" --conda_base "mamba" --print_only
```
`--n_jobs` - amount of parallel jobs and requested cpu cores if run using sbatch

`--total_memory_gb` - total RAM memory requested, either per HPC node or locally

`--qos` - QOS option passed to sbatch script if run on HPC. This is HPC dependent

`--conda_base` - type of conda installation if using sbatch
  Accepted values are "mamba", "mambaforge", "conda", "miniconda", "anaconda", "miniforge", "miniforge3"
  If loaded as a module on HPC specify "module <module_name>", e.g. "module mamba"
  You can also provide custom path to your conda installation e.g. "/custom/path/to/conda.sh"

`--print_only` - writes snakemake, sbatch and qsub scripts into your_cell_level_directory/snakemake/

### Run mapping from generated script files
```shell
bash your_cell_level_directory/snakemake/qsub/snakemake_cmd.txt
```
```shell
bash your_cell_level_directory/snakemake/sbatch/sbatch-serial-qos.sh
```

### Run mapping interactively in the current shell
```shell
yap mapping --output_dir "your_cell_level_directory" \
    --config_path "m3c_config.ini" --n_jobs 14 --total_memory_gb 32
```

## Map already-demultiplexed cell-level FASTQ files from external sources

Use `mapping-cell-fastq` when you have already-demultiplexed cell-level FASTQ files that are not organized in pipeline expected directory structure, e.g. from external sources. Cells are randomly shuffled into groups and symlinked into the expected directory structure, then the mapping pipeline is prepared for each group.

FASTQ files must follow the canonical naming format `{plate}-{well}-R1.fq.gz` / `R2.fq.gz`. If your files come from an external source with a different naming convention, run `link-cell-fastq` first (see above).

```shell
yap mapping-cell-fastq \
    --output_dir "your_cell_level_directory" \
    --config_path "m3c_config.ini" \
    --fastq_pattern "cell_fastq/*.fq.gz" \
    --cells_per_group 64 \
    --n_jobs 64 \
    --total_memory_gb 128 \
    --qos "serial" \
    --conda_base "mamba"
```

**Parameters:**
- `-o` / `--output_dir`: (**required**) Output directory. Must not already exist.
- `-config` / `--config_path`: (**required**) Path to the mapping config `.ini` file (see `yap default-mapping-config`).
- `-fq` / `--fastq_pattern`: (**required**) Glob pattern matching all cell-level FASTQ files. **Must be quoted** to prevent shell expansion, e.g. `"cell_fastq/*.fq.gz"`.
- `--cells_per_group`: Number of cells per batch subdirectory. Cells are randomly shuffled into groups of this size. Default: 64.
- `--n_jobs`: Number of parallel jobs per group.
- `--total_memory_gb`: Total RAM available.
- `--qos`: QOS for sbatch script.
- `--conda_base`: Conda installation type. See `--conda_base` options under **Mapping** above.

This will create mapping scripts, that you can execute either locally or by submitting HPC jobs:

```shell
bash your_cell_level_directory/snakemake/qsub/snakemake_cmd.txt
```
```shell
bash your_cell_level_directory/snakemake/sbatch/sbatch-serial-qos.sh
```


# Background on library preparation and structure

## Library preparation
384 random index primers are used, so each well on 384 well plates receives separate unique random index

<img src="doc/molecularsteps.png" width="600px">

## Version 2 barcoding strategy
Please note 1 single PCR index primer can be used for the whole 384 well plate, in this case there will be only one multiplex group
<img src="doc/v2barcode.png" width="800px">


## Library structure
<img src="doc/library.png" width="300px">


# Workflow of m3c run
<img src="doc/snm3c_dag.svg" title="DAG for snm3c" width="800px">

# QC metrics

| Field | Purpose |
|-------|---------|
| **cell** | The unique identifier for each cell |
| **InputReadPairs** | Total number of reads input into Cutadapt for processing |
| **InputReadPairsBP** | Total number of base pairs in the input reads |
| **TrimmedReadPairs** | Number of reads that were written to the output after trimming and filtering |
| **R1WithAdapters** | Number of reads that contained adapter sequences that were detected and trimmed |
| **R1QualTrimBP** | Number of R1 base pairs trimmed due to low base quality |
| **R1TrimmedReadsBP** | Number of R1 base pairs remaining after adapter and quality trimming |
| **R2WithAdapters** | Number of reads from the second mate in paired-end reads that contained adapter sequences and were trimmed |
| **R2QualTrimBP** | Number of R2 base pairs trimmed due to low base quality |
| **R2TrimmedReadsBP** | Number of R2 base pairs remaining after adapter and quality trimming |
| **ReadsUniqueMapped** | **`mc`** specific. Number of uniquely mapped reads, calculated from Hisat3n report: 2 * PEUniqueMappedReadPairs (Read Pairs Aligned concordantly 1 time) + 2 * PEDiscordantlyUniqueMappedReadPairs (Read Pairs Aligned discordantly 1 time) + SEUniqueMappedReads (Unpaired Reads Aligned 1 time) |
| **ReadsUniqueMappingRate** | **`mc`** specific. Rate of unique read mapping: (ReadsUniqueMapped / (ReadPairsMappedInPE (Total Read Pairs) * 2) * 100) |
| **ReadsMultiMapped** | **`mc`** specific. Number of multimapped reads:  (2 * PEMultiMappedReadPairs (Read Pairs Aligned concordantly >1 times)) + SEMultiMappedReads (Unpaired Reads Aligned >1 times) |
| **ReadsMultiMappingRate** | **`mc`** specific. Rate of multimapping: (ReadsMultiMapped / (ReadPairsMappedInPE (Total Read Pairs) * 2) * 100) |
| **ReadsOverallMappingRate** | **`mc`** specific. Rate of mapping for all reads (ReadsUniqueMappingRate + ReadsMultiMappingRate) |
| **DNAReadsUniqueMapped** | **`mct`** specific. Number of uniquely mapped DNA reads, calculated from Hisat3n DNA alignment report: 2 * PEUniqueMappedReadPairs + 2 * PEDiscordantlyUniqueMappedReadPairs + SEUniqueMappedReads |
| **DNAReadsUniqueMappingRate** | **`mct`** specific. Rate of unique DNA read mapping: (DNAReadsUniqueMapped / (ReadPairsMappedInPE * 2) * 100) |
| **DNAReadsMultiMapped** | **`mct`** specific. Number of multimapped DNA reads: (2 * PEMultiMappedReadPairs) + SEMultiMappedReads |
| **DNAReadsMultiMappingRate** | **`mct`** specific. Rate of DNA multimapping: (DNAReadsMultiMapped / (ReadPairsMappedInPE * 2) * 100) |
| **DNAReadsOverallMappingRate** | **`mct`** specific. Overall DNA mapping rate (DNAReadsUniqueMappingRate + DNAReadsMultiMappingRate) |
| **RNAReadsUniqueMapped** | **`mct`** specific. Number of uniquely mapped RNA reads from Hisat3n RNA alignment report: 2 * PEUniqueMappedReadPairs + 2 * PEDiscordantlyUniqueMappedReadPairs + SEUniqueMappedReads |
| **RNAReadsUniqueMappingRate** | **`mct`** specific. Rate of unique RNA read mapping: (RNAReadsUniqueMapped / (ReadPairsMappedInPE * 2) * 100) |
| **RNAReadsMultiMapped** | **`mct`** specific. Number of multimapped RNA reads: (2 * PEMultiMappedReadPairs) + SEMultiMappedReads |
| **RNAReadsMultiMappingRate** | **`mct`** specific. Rate of RNA multimapping: (RNAReadsMultiMapped / (ReadPairsMappedInPE * 2) * 100) |
| **RNAReadsOverallMappingRate** | **`mct`** specific. Overall RNA mapping rate (RNAReadsUniqueMappingRate + RNAReadsMultiMappingRate) |
| **WholeReadsUniqueMapped** | **`m3C`** specific. Same calculation as in **ReadsUniqueMapped**. Number of uniquely mapped reads from initial hisat3n alignment, before m3c-specific restriction-site splitting. Each unmapped read from this step is split into fragments and remapped (see R1/R2SplitReads). Underreports true mapping rate since it misses split-read mappings |
| **WholeReadsMultiMapped** | **`m3C`** specific. Same calculation as in **ReadsMultiMapped**. Number of multimapped reads from initial hisat3n alignment, before restriction-site splitting |
| **WholeReadsUniqueMappingRate** | **`m3C`** specific. Same calculation as in **ReadsUniqueMappingRate**. Rate of unique mapping from initial hisat3n alignment. Underreports true unique mapping rate since it misses split-read mappings. See UniqueClusterMappingRate for the correct overall rate |
| **WholeReadsMultiMappingRate** | **`m3C`** specific. Same calculation as in **ReadsMultiMappingRate**. Rate of multimapping from initial hisat3n alignment |
| **WholeReadsOverallMappingRate** | **`m3C`** specific. Same calculation as in **ReadsOverallMappingRate**. Overall mapping rate from initial hisat3n alignment. Underreports true overall mapping rate since it misses split-read mappings |
| **R1SplitReadsUniqueMapped** | **`m3C`** specific. Number of uniquely mapped R1 reads that failed initial hisat3n paired-end mapping, were split at restriction cut sites, and remapped as single-end reads |
| **R1SplitReadsUniqueMappingRate** | **`m3C`** specific. Rate of unique mapping for R1 split reads |
| **R1SplitReadsMultiMapped** | **`m3C`** specific. Number of multimapped R1 split reads |
| **R1SplitReadsMultiMappingRate** | **`m3C`** specific. Rate of multimapping for R1 split reads |
| **R1SplitReadsOverallMappingRate** | **`m3C`** specific. Overall mapping rate for R1 split reads |
| **R2SplitReadsUniqueMapped** | **`m3C`** specific. Number of uniquely mapped R2 reads that failed initial hisat3n paired-end mapping, were split at restriction cut sites, and remapped as single-end reads |
| **R2SplitReadsUniqueMappingRate** | **`m3C`** specific. Rate of unique mapping for R2 split reads |
| **R2SplitReadsMultiMapped** | **`m3C`** specific. Number of multimapped R2 split reads |
| **R2SplitReadsMultiMappingRate** | **`m3C`** specific. Rate of multimapping for R2 split reads |
| **R2SplitReadsOverallMappingRate** | **`m3C`** specific. Overall mapping rate for R2 split reads |
| **UniqueMappedClusters** | **`m3C`** specific.  Number of unique read clusters (QNAMEs) in the BAM file, including both initial reads and cut-site split reads. Correctly reports total mapped clusters, in contrast to WholeReadsUniqueMapped missing m3c split-reads and R1/R2SplitReadsUniqueMapped missing initial whole reads |
| **UniqueMappedR1** | **`m3C`** specific. Number of unique mapped R1 in the BAM file, including both initial reads and cut-site split reads. Correctly reports total mapped R1, in contrast to WholeReadsUniqueMapped missing m3c R1 split-reads and R1SplitReadsUniqueMapped missing initial whole R1 reads |
| **UniqueMappedR2** | **`m3C`** specific. Number of unique mapped R2 in the BAM file, including both initial reads and cut-site split reads. Correctly reports total mapped R2, in contrast to WholeReadsUniqueMapped missing m3c R2 split-reads and R2SplitReadsUniqueMapped missing initial whole R2 reads |
| **UniqueClusterMappingRate** | **`m3C`** specific. UniqueMappedClusters / TrimmedReadPairs × 100. BAM-derived overall unique mapping rate. Correctly reports cluster mapping rate during m3c |
| **UniqueMappedR1Rate** | **`m3C`** specific. UniqueMappedR1 / TrimmedReadPairs × 100. BAM-derived R1 unique mapping rate. Correctly reports R1 mapping rate during m3c |
| **UniqueMappedR2Rate** | **`m3C`** specific. UniqueMappedR2 / TrimmedReadPairs × 100. BAM-derived R2 unique mapping rate. Correctly reports R2 mapping rate during m3c |
| **CisContacts** | **`m3C`** specific. R1 and R2 from the same read pair map to different loci on the same chromosome, ≥ `min_gap` bp apart. These are the primary source of cis chromatin contacts |
| **CisCutContacts** | **`m3C`** specific. Two adjacent fragments from the same read (either R1 or R2) split at a restriction enzyme cut site (NlaIII/DpnII) both map to the same chromosome ≥ `min_gap` bp apart. Cut-site-confirmed cis contacts — the restriction site is directly sequenced within the read. Less frequent than CisContacts since read must physically span a cut site |
| **CisMultiContacts** | **`m3C`** specific.  A R1-R2 read pair producing 3 or more distinct mapping locations (≥ `min_gap` bp apart or on different chromosomes), where at least one read was split at a restriction site. This specific R1-R2 contact is between two locations on the same chromosome. Other contacts from the same read pair are classified independently and may either be CisCutMultiContacts, TransMultiContacts, or TransCutMultiContacts |
| **CisCutMultiContacts** | **`m3C`** specific. Adjacent restriction cut-site-split fragments from the same read, part of a R1-R2 read pair or single read producing 3 or more distinct mapping locations (≥ `min_gap` bp apart or on different chromosomes). This specific contact is between two adjacent fragments on the same chromosome. Co-occurs with CisMultiContacts, TransMultiContacts, or TransCutMultiContacts from the same read pair depending on where the other locations land |
| **TransContacts** | **`m3C`** specific. R1 and R2 from the same read pair map to different chromosomes. Inter-chromosomal contacts |
| **TransCutContacts** | **`m3C`** specific. Two adjacent fragments from the same read split at a restriction enzyme cut site map to different chromosomes. Cut-site-confirmed trans contacts |
| **TransMultiContacts** | **`m3C`** specific. A R1-R2 read pair producing 3 or more distinct mapping locations where this specific R1-R2 contact is between two locations on different chromosomes. Requires at least one read to be split at a restriction site. The same read pair simultaneously produces other contacts classified independently as CisMultiContacts, CisCutMultiContacts, or TransCutMultiContacts depending on chromosome and cut-site adjacency of each pair |
| **TransCutMultiContacts** | **`m3C`** specific. Adjacent restriction cut-site-split fragments from the same read, that is a member of R1-R2 read pair or single read mapping to 3 or more separate locations on different chromosomes. Co-occurs with TransMultiContacts, CisMultiContacts, or CisCutMultiContacts from the same read pair depending on where the other locations land |
| **ChimericContacts** | **`m3C`** specific. Two fragments from the same read that are not adjacent at a restriction cut site. Indicates a chimeric DNA molecule — two unrelated sequences ligated without restriction digestion during library prep. Tracked in QC summary, but excluded from final contact files |
| **NoContacts** | **`m3C`** specific. Read pairs where all mapping locations fell within `min_gap` bp of each other on the same chromosome — no contact could be called. Tracked in QC summary, but excluded from final contact files |
| **MappedFragments** | **`m3C`** specific. Total number of read pair names (fragments) processed by contact calling, regardless of outcome (contact, no-contact, chimeric). Used as denominator for CisContactsRatio, TransContactsRatio, MultiContactsRatio |
| **DeduppedContacts** | **`m3C`** specific. Total number of deduplicated contacts |
| **ContactsDeduplicationRate** | **`m3C`** specific. (input_contacts - dedup_contacts) / (input_contacts + 0.00001) |
| **TotalCisContacts** | **`m3C`** specific. Total number of cis contacts |
| **TotalTransContacts** | **`m3C`** specific. Total number of trans contacts |
| **TotalMultiContacts** | **`m3C`** specific. Total number of multi contacts |
| **CisContactsRatio** | **`m3C`** specific. Ratio of cis contacts: TotalCisContacts / MappedFragments |
| **TransContactsRatio** | **`m3C`** specific. Ratio of trans contacts: TotalTransContacts / MappedFragments |
| **MultiContactsRatio** | **`m3C`** specific. Ratio of multi contacts: TotalMultiContacts / MappedFragments |
| **UniqueAlignFinalReads** | **`mc/m3C`** specific. Final unique mapped total reads after picard deduplication |
| **UniqueAlignDuplicatedReads** | **`mc/m3C`** specific. Paired and unpaired duplicated reads from unique alignment |
| **UniqueAlignPCRDuplicationRate** |  **`mc/m3C`** specific. The percentage of reads identified as PCR duplicates: (1 - UniqueAlignFinalReads / (UniqueAlignFinalReads + UniqueAlignDuplicatedReads)) * 100 |
| **MultiAlignFinalReads** | **`mc/m3C`** specific. Final multimapped total reads after picard deduplication |
| **MultiAlignDuplicatedReads** | **`mc/m3C`** specific. Paired and unpaired duplicated reads from multi alignment |
| **MultiAlignPCRDuplicationRate** | **`mc/m3C`** specific. The percentage of multi-mapped reads identified as PCR duplicates: (1 - MultiAlignFinalReads / (MultiAlignFinalReads + MultiAlignDuplicatedReads)) * 100 |
| **DNAUniqueAlignFinalReads** | **`mct`** specific. Final unique mapped DNA reads after picard deduplication |
| **DNAUniqueAlignDuplicatedReads** | **`mct`** specific. Paired and unpaired duplicated reads from unique DNA alignment |
| **DNAUniqueAlignPCRDuplicationRate** | **`mct`** specific. The percentage of unique DNA reads identified as PCR duplicates: (1 - DNAUniqueAlignFinalReads / (DNAUniqueAlignFinalReads + DNAUniqueAlignDuplicatedReads)) * 100 |
| **DNAMultiAlignFinalReads** | **`mct`** specific. Final multimapped DNA reads after picard deduplication |
| **DNAMultiAlignDuplicatedReads** | **`mct`** specific. Paired and unpaired duplicated reads from multi DNA alignment |
| **DNAMultiAlignPCRDuplicationRate** | **`mct`** specific. The percentage of multi DNA reads identified as PCR duplicates: (1 - DNAMultiAlignFinalReads / (DNAMultiAlignFinalReads + DNAMultiAlignDuplicatedReads)) * 100 |
| **UniqueAlignFinalDNAReads** | **`mct`** specific. Number of uniquely mapped DNA reads selected from unique alignment BAM based on mCH fraction threshold |
| **UniqueAlignSelectedDNAReadsRatio** | **`mct`** specific. Fraction of reads passing DNA selection from unique alignment BAM |
| **MultiAlignFinalDNAReads** | **`mct`** specific. Number of multimapped DNA reads selected from multi alignment BAM based on mCH fraction threshold |
| **MultiAlignSelectedDNAReadsRatio** | **`mct`** specific. Fraction of reads passing DNA selection from multi alignment BAM |
| **FinalRNAReads** | **`mct`** specific. Number of RNA reads selected from RNA alignment BAM based on mCH fraction threshold |
| **SelectedRNAReadsRatio** | **`mct`** specific. Fraction of reads passing RNA selection from RNA alignment BAM |
| **AssignedRNAReads** | **`mct`** specific. Number of reads assigned to genomic features (genes) by featureCounts |
| **UnassignedRNAReads** | **`mct`** specific. Total number of reads not assigned to any genomic feature |
| **AssignedRNAReadsRate** | **`mct`** specific. AssignedRNAReads / (AssignedRNAReads + UnassignedRNAReads) |
| **mCCC** | Total methylated cytosines in the CCC context |
| **mCG** | Total methylated cytosines in the CG context |
| **mCH** | Total methylated cytosines in the CH context |
| **CCCCov** | Total covered cytosines in the CCC context |
| **CGCov** | Total covered cytosines in the CG context |
| **CHCov** | Total covered cytosines in the CH context |
| **mCCCFrac** | Fraction of methylated cytosines (mCCC) divided by covered cytosines (CCCCov) in the CCC context |
| **mCGFrac** | Fraction of methylated cytosines (mCG) divided by covered cytosines (CGCov) in the CG context |
| **mCHFrac** | Fraction of methylated cytosines (mCH) divided by covered cytosines (CHCov) in the CH context |
| **GenomeBreadth** | Fraction of all genome nucleotide positions covered by at least one read (breadth of coverage) |
| **LambdamCGFrac** | Fraction of methylated CG cytosines divided by covered CG cytosines on lambda phage spike-in DNA (chrL). Should be ~1 if enzymatic CG methylation was complete. Drop below 1 indicates over-conversion when mCG are stripped of methyl group, or insufficient enzymatic methylation of lambda DNA |
| **LambdaCGCov** | Total covered CG cytosines on lambda spike-in DNA (chrL) |
| **LambdamCHFrac** | Fraction of methylated CH (CA + CT + CC) cytosines divided by covered CH cytosines on lambda spike-in DNA (chrL). Should be ~0 since lambda CH in unmethylated, high values indicate underconversion of unmethylated CH into TH |
| **LambdaCHCov** | Total covered CH cytosines on lambda spike-in DNA (chrL) |
