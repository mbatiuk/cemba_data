import pathlib
import re
from collections import defaultdict
import pandas as pd
import yaml
from pysam import TabixFile
from ..utilities import get_configuration


def _read_yaml_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def _read_ini_config(config_path):
    return get_configuration(config_path)


def read_mapping_config(cwd: str = '.'):
    tried = []
    yaml_path = None
    for name in ['config', 'mapping_config']:
        for config_dir in [cwd, f'{cwd}/..']:
            for suffix in ['yaml', 'yml']:
                path = f'{config_dir}/{name}.{suffix}'
                tried.append(path)
                if pathlib.Path(path).exists():
                    yaml_path = path
    default_path = '~/mapping_config.yaml'
    if pathlib.Path(default_path).exists():
        yaml_path = default_path

    ini_path = None
    for name in ['config', 'mapping_config']:
        for config_dir in [cwd, f'{cwd}/..']:
            path = f'{config_dir}/{name}.ini'
            tried.append(path)
            if pathlib.Path(path).exists():
                ini_path = path

    if yaml_path is not None:
        config = _read_yaml_config(yaml_path)
    elif ini_path is not None:
        config = _read_ini_config(ini_path)
    else:
        config = {}
    return config


def get_allc_lambda_frac(allc_list, num_upstr_bases):
    """Compute bisulfite conversion rate from lambda spike-in DNA (chrL)."""
    num_upstr_bases = int(num_upstr_bases)
    records = {}
    for path in allc_list:
        mc_counts = defaultdict(int)
        cov_counts = defaultdict(int)
        with TabixFile(str(path)) as allc:
            cell = pathlib.Path(path).name.split('.')[0]
            try:
                for line in allc.fetch('chrL'):
                    chrom, pos, strand, context, mc, cov, _ = line.split('\t')
                    context = context[num_upstr_bases:num_upstr_bases + 2]
                    mc_counts[context] += int(mc)
                    cov_counts[context] += int(cov)
                df = pd.DataFrame({'mc': pd.Series(mc_counts), 'cov': pd.Series(cov_counts)})
                df = df.reindex(['CG', 'CC', 'CT', 'CA']).fillna(0)
                cy_cov = df.loc['CT', 'cov'] + df.loc['CC', 'cov']
                cy_frac = (
                    (df.loc['CT', 'mc'] + df.loc['CC', 'mc']) / cy_cov
                    if cy_cov > 0 else 0
                )
                records[cell] = {'LambdaCYFrac': cy_frac, 'LambdaCYCov': cy_cov}
            except ValueError:
                records[cell] = {'LambdaCYFrac': 0, 'LambdaCYCov': 0}
    return pd.DataFrame(records).T
