import os,sys
import pandas as pd
import pathlib
import re
import glob
from cemba_data.demultiplex import _parse_index_fasta,_read_cutadapt_result, \
get_fastq_info, get_random_index

MODULE_DIR = pathlib.Path(workflow.snakefile).parent.absolute()

default_config={
    'total_read_pairs_min':1,
    'total_read_pairs_max':10000000
}

for key in default_config:
    if key not in config:
        config[key]=default_config[key]

fq_dir = pathlib.Path(config["fq_dir"]).absolute()

outdir=config.get("outdir","mapping")

cells_per_group = int(config.get("cells_per_group", 64))

df_fq=get_fastq_info(fq_dir,outdir)
uid_fastqs_dict=df_fq.loc[:,['uid','R1','R2']].groupby('uid').agg(lambda x:x.tolist()).to_dict(orient='index')

df_index=get_random_index(df_fq.uid.unique().tolist(), cells_per_group, outdir)

rule summary_demultiplex:
    input:
        stats=expand(outdir+"/{uid}_fq_unknown/demultiplex.stats.txt",
                                    uid=df_fq.uid.unique().tolist())
    output:
        csv=outdir+"/stats/demultiplex.stats.csv",
        fq_info=outdir+"/stats/fastq_info.tsv",
        index_info=outdir+"/stats/index_info.tsv"
    params:
        stat_dir=os.path.join(outdir,"stats")
    run:
        shell(f"mkdir -p {params.stat_dir}")
        df_fq.to_csv(output.fq_info,sep='\t',index=False)
        df_index.to_csv(output.index_info,sep='\t',index=False)
        random_index_fasta_path = MODULE_DIR / 'random_index_v2.fa'
        index_seq_dict = _parse_index_fasta(random_index_fasta_path)
        index_name_dict = {v: k for k, v in index_seq_dict.items()}
        stat_list = []
        for path in glob.glob(outdir + "/*_fq_unknown/demultiplex.stats.txt"):
            uid = path.split('/')[-2].replace('_fq_unknown', '')
            single_df=_read_cutadapt_result(path)
            single_df['uid'] = uid
            single_df['index_name'] = single_df['Sequence'].map(index_name_dict)
            assert single_df['index_name'].isna().sum() == 0
            stat_list.append(single_df)
        df_stats = pd.concat(stat_list)
        # look up group from actual fastq dirs created by run_demultiplex
        def get_group_from_dirs(uid, index_name):
            for d in glob.glob(outdir + f"/{uid}Group*/fastq/{uid}-{index_name}-R1.fq.gz"):
                return d.split('/')[-3].replace(uid + 'Group', '')
            return 'NA'
        df_stats['group'] = df_stats.apply(
            lambda r: get_group_from_dirs(r['uid'], r['index_name']), axis=1)
        df_stats = df_stats.loc[df_stats['group'] != 'NA']
        df_stats['cell_id'] = df_stats['uid'] + '-' + df_stats['index_name']
        df_cell=df_stats.groupby('cell_id').agg({
                                        'Trimmed':'sum',
                                        'TotalPair':'sum',
                                        'index_name':lambda i: i.unique()[0],
                                        'uid':lambda i: i.unique()[0]})
        df_cell.rename(columns={'Trimmed': 'CellInputReadPairs',
                                'TotalPair': 'MultiplexedTotalReadPairs',
                                'index_name': 'IndexName',
                                'uid': 'UID'},inplace=True)
        df_cell['CellBarcodeRate'] = df_cell['CellInputReadPairs'] / df_cell['MultiplexedTotalReadPairs']
        df_cell.to_csv(output.csv)
        if os.path.exists(os.path.join(outdir,"fastq_info.txt")):
            os.remove(os.path.join(outdir,"fastq_info.txt"))
        if os.path.exists(os.path.join(outdir,"random_index.txt")):
            os.remove(os.path.join(outdir,"random_index.txt"))
        # move cutadapt logs to stats/cutadapt/
        shell(f"mkdir -p {params.stat_dir}/cutadapt")
        for path in glob.glob(outdir + "/*_fq_unknown/demultiplex.stats.txt"):
            uid = path.split('/')[-2].replace('_fq_unknown', '')
            os.rename(path, os.path.join(outdir, "stats", "cutadapt", uid + ".cutadapt.stats.txt"))

rule merge_lanes:
    input:
        fqs=lambda wildcards: uid_fastqs_dict[wildcards.uid][wildcards.read_type]

    output:
        fq=temp(outdir+"/{uid}_fq_unknown/{read_type}.fq.gz")

    params:
        outdir=lambda wildcards: outdir+f"/{wildcards.uid}_fq_unknown"

    run:
        shell(f"mkdir -p {params.outdir}")
        shell(f"cat {input.fqs} > {output.fq}")

rule run_demultiplex:
    input:
        R1 = lambda wildcards: outdir+f"/{wildcards.uid}_fq_unknown/R1.fq.gz",
        R2 = lambda wildcards: outdir+f"/{wildcards.uid}_fq_unknown/R2.fq.gz"

    output:
        stats_out = outdir+"/{uid}_fq_unknown/demultiplex.stats.txt",

    params:
        random_index_fa = MODULE_DIR / 'random_index_v2.fa',
        outdir=lambda wildcards: outdir+f"/{wildcards.uid}_fq_unknown/demultiplex",
        R1=lambda wildcards: outdir+f"/{wildcards.uid}_fq_unknown/demultiplex/{'{{name}}'}-R1.fq.gz",
        R2=lambda wildcards: outdir+f"/{wildcards.uid}_fq_unknown/demultiplex/{'{{name}}'}-R2.fq.gz"
    run:
        shell(f"mkdir -p {params.outdir}")
        shell(f"cutadapt -e 0.13 --no-indels -g file:{params.random_index_fa} -o  {params.R1} -p {params.R2} {input.R1} {input.R2} > {output.stats_out}")
         # for the reads startswith random index present in random_index_fa, will be taken and write into 1 fastq (1 cell),
         # cut the left 8 bp sequence and add the random index name (A2, P24) into the cell fastq name.
         # one uid will be broken down into 384 cells.

        # remove temporary fastq files
        os.remove(input.R1)
        os.remove(input.R2)

        # parse demultiplex.stats.txt for cell_qc
        random_index_fasta_path = MODULE_DIR / 'random_index_v2.fa'
        index_seq_dict = _parse_index_fasta(random_index_fasta_path)
        index_name_dict = {v: k for k, v in index_seq_dict.items()}
        single_df = _read_cutadapt_result(output.stats_out)
        single_df['index_name'] = single_df['Sequence'].map(index_name_dict)
        removed_index_names = single_df.loc[
            (single_df['Trimmed'] < int(config['total_read_pairs_min'])) |
            (single_df['Trimmed'] > int(config['total_read_pairs_max']))
        ].index_name.unique().tolist()

        # collect valid R1 fastqs, sort for consistent group assignment
        all_fqs = glob.glob(outdir+f"/{wildcards.uid}_fq_unknown/demultiplex/*-R1.fq.gz")
        valid_r1 = sorted([f for f in all_fqs
                           if 'unknow' not in os.path.basename(f).lower()
                           and os.path.basename(f).split('-')[0] not in removed_index_names])

        # move demultiplexed fastq files to {plate}Group{n}/fastq/
        for i, r1_fq in enumerate(valid_r1):
            group = (i // cells_per_group) + 1
            index_name = os.path.basename(r1_fq).split('-')[0]
            new_uid = wildcards.uid + 'Group' + str(group)
            for read_type in ['R1', 'R2']:
                src = outdir+f"/{wildcards.uid}_fq_unknown/demultiplex/{index_name}-{read_type}.fq.gz"
                dst = outdir+f"/{new_uid}/fastq/{wildcards.uid}-{index_name}-{read_type}.fq.gz"
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                os.rename(src, dst)

        # move unknown fastqs to parent _fq_unknown dir
        for fq in glob.glob(outdir+f"/{wildcards.uid}_fq_unknown/demultiplex/unknown-R*.fq.gz"):
            os.rename(fq, outdir+f"/{wildcards.uid}_fq_unknown/{os.path.basename(fq)}")

        # move filtered-out cell fastqs (failed min/max thresholds) to filtered/
        remaining = glob.glob(outdir+f"/{wildcards.uid}_fq_unknown/demultiplex/*-R*.fq.gz")
        if remaining:
            filtered_dir = outdir+f"/{wildcards.uid}_fq_unknown/filtered"
            os.makedirs(filtered_dir, exist_ok=True)
            for fq in remaining:
                os.rename(fq, os.path.join(filtered_dir, os.path.basename(fq)))

        # remove now-empty demultiplex subdir
        demux_dir = outdir+f"/{wildcards.uid}_fq_unknown/demultiplex"
        if os.path.isdir(demux_dir) and not os.listdir(demux_dir):
            os.rmdir(demux_dir)