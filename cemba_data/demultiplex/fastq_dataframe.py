"""
Generate raw FASTQ dataframe based on fixed name pattern
"""
import os
import glob
import logging
import pathlib
import pandas as pd

# logger
log = logging.getLogger()


def _parse_fastq_path(path):
    """
    UID pattern: {sample_id_prefix}-{plate}-{multiplex_group}-{barcode_name}
    FASTQ name pattern:
    {sample_id_prefix}-{plate}-{multiplex_group}-{barcode_name}_{internal_info}_{lane}_{read_type}_{internal_info}.fastq.gz
    """
    try:
        *_, plate, multiplex_group, multi_field = os.path.basename(path).split('-')
        primer_name, *_, lane, read_type, _ = multi_field.split('_')
        try:
            assert primer_name[0] in 'ABCDEFGHIJKLMNOP'
            assert int(primer_name[1:]) in list(range(1, 25))
            assert int(multiplex_group) in list(range(1, 7))
            assert lane in {'L001', 'L002', 'L003', 'L004', 'L005', 'L006', 'L007', 'L008'}
        except AssertionError:
            raise ValueError
    except ValueError:
        raise ValueError(f'Found unknown name pattern in path {path}')
    name_dict = dict(plate=plate,
                     multiplex_group=multiplex_group,
                     primer_name=primer_name,
                     lane=lane,
                     read_type=read_type,
                     fastq_path=path,
                     uid=f'{plate}-{multiplex_group}-{primer_name}')
    return pd.Series(name_dict)

def make_fastq_dataframe(file_path, output_path=None):
    """
    Generate fastq_dataframe for pipeline input.

    Parameters
    ----------
    file_path
        Accept 1. path pattern containing wildcard, 2. path list,
        3. path of one file containing all the paths.
    output_path
        output path of the fastq dataframe

    Returns
    -------
    fastq_dataframe for pipeline input.
    """
    if isinstance(file_path, str) and ('*' in file_path):
        file_path = [str(pathlib.Path(p).absolute()) for p in glob.glob(file_path)]
    elif isinstance(file_path, list):
        pass
    else:
        with open(file_path) as f:
            file_path = [line.strip() for line in f]
    log.info(f'{len(file_path)} FASTQ file paths in input')

    fastq_data = [_parse_fastq_path(path) for path in file_path]
    fastq_df = pd.DataFrame(fastq_data)
    log.info(f'{fastq_df.shape[0]} valid fastq names.')
    if fastq_df.shape[0] == 0:
        log.info('No fastq name remained, check if the name pattern is correct.')
        return None

    # make sure UID is unique
    for _, df in fastq_df.groupby(['lane', 'read_type']):
        if df['uid'].unique().size != df['uid'].size:
            raise ValueError('UID column is not unique.')
    if output_path is not None:
        fastq_df.to_csv(output_path, index=False)
    return fastq_df
