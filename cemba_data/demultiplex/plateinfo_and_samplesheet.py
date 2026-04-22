"""
Contain codes about parse plate info and generate sample sheet
"""

import pathlib
import re
from collections import OrderedDict

import pandas as pd

import cemba_data

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])

# the Illumina sample sheet header used by Ecker Lab
with open(PACKAGE_DIR / 'files/sample_sheet_header.txt') as _f:
    SAMPLESHEET_DEFAULT_HEADER = _f.read()

SECTIONS = ['[CriticalInfo]', '[LibraryInfo]', '[PlateInfo]']

LIMITED_CHOICES = {
    'input_plate_size': [384, '384'],
    'primer_quarter': ['Set1_Q1', 'Set1_Q2', 'Set1_Q3', 'Set1_Q4',
                       'SetB_Q1', 'SetB_Q2', 'SetB_Q3', 'SetB_Q4']}

CRITICAL_INFO_KEYS = ['input_plate_size', 'pool_id', 'tube_label', 'email']

BARCODE_TABLE_PATH = PACKAGE_DIR / 'files/V2_i7_i5_index.tsv'


def _clean_str_for_path(str_in):
    # replace special char with _
    str_out = re.sub('[^a-zA-Z0-9]', '_', str_in.strip())
    return str_out


def _get_kv_pair(line):
    try:
        k, v = line.split('=')
        if k == 'email':
            return k, v
        else:
            return _clean_str_for_path(k), _clean_str_for_path(v)
    except ValueError:
        raise ValueError(f'Each key=value line must contain a "=" to separate key and value. Got {line}')


def _read_plate_info(plate_info_path):
    """Parse the plate info file"""
    cur_section = ''
    cur_section_id = -1

    critical_info = {}
    library_info = OrderedDict()
    plate_header = True
    plate_info = []

    with open(plate_info_path) as f:
        for line in f:
            line = line.strip('\n')
            if line == '' or line.startswith('#'):
                continue

            if line.startswith('['):
                cur_section_id += 1
                if line == SECTIONS[cur_section_id]:
                    cur_section = line
                else:
                    raise ValueError(
                        f'Section name and order must be [CriticalInfo] [LibraryInfo] [PlateInfo], '
                        f'got {line} at No.{cur_section_id + 1} section.')
            elif cur_section == '[CriticalInfo]':
                k, v = _get_kv_pair(line)
                if k not in CRITICAL_INFO_KEYS:
                    raise ValueError(f'Unknown key {k} in [CriticalInfo]')
                else:
                    critical_info[k] = v
            elif cur_section == '[LibraryInfo]':
                k, v = _get_kv_pair(line)
                if (k in critical_info.keys()) or (k in library_info.keys()):
                    raise ValueError(f'Found duplicated key {k}')
                else:
                    library_info[k] = v
            elif cur_section == '[PlateInfo]':
                ll = re.split(r' |,|\t', line)
                if plate_header:
                    plate_header = False
                plate_info.append(ll)
            else:
                raise ValueError(f'Got a malformed line: {line}')

    for k in CRITICAL_INFO_KEYS:
        if k not in critical_info:
            raise ValueError(f'[CriticalInfo] missing key-value pair "{k}"')

    header = plate_info[0]
    plate_info = pd.DataFrame(plate_info[1:], columns=plate_info[0])
    for k, v in library_info.items():
        if k in plate_info.columns:
            raise ValueError(f'Found duplicated key {k} between [PlateInfo] and [LibraryInfo]')
        plate_info[k] = v

    n_plate_info_fix_col = 3
    if plate_info.set_index(['plate_id', 'multiplex_group']).index.duplicated().sum() != 0:
        raise ValueError('Found duplicated plate_id, multiplex_group combination in [PlateInfo] section.')

    col_order = header[:n_plate_info_fix_col] + list(library_info.keys()) + header[n_plate_info_fix_col:]
    plate_info = plate_info[col_order].copy()
    plate_info['sample_id_prefix'] = plate_info.apply(
        lambda i: '-'.join(i[n_plate_info_fix_col:].astype(str).tolist()), axis=1)

    for k, v in critical_info.items():
        if k in plate_info.columns:
            raise ValueError(f'Found duplicated key {k}')
        plate_info[k] = v

    plate_info = plate_info.loc[~plate_info['primer_name'].duplicated(), :].copy()
    return critical_info, plate_info


def reverse_comp(sequence):
    rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
               'a': 't', 't': 'a', 'c': 'G', 'g': 'C',
               'n': 'n', 'N': 'N'}
    _seq = ''.join([rc_dict[s] for s in sequence[::-1]])
    return _seq


def _plate_384_random_index_384(plate_info, barcode_table, i5_reverse_comp=False):
    """
    UID pattern: {sample_id_prefix}-{plate}-{multiplex_group}-{barcode_name}

    If i5_reverse_comp, use the reverse complement i5 sequence,
    needed for NovaSeq v1.5 S4 kit.
    """
    records = []

    for _, row in plate_info.iterrows():
        plate = _clean_str_for_path(row['plate_id'])

        barcode_name = row['primer_name']
        cur_row = barcode_table.loc[barcode_name]
        i5_barcode = cur_row['i5_index_sequence']
        i7_barcode = cur_row['i7_index_sequence']
        sample_id_prefix = row['sample_id_prefix']
        multiplex_group = row['multiplex_group']
        sample_id = f'{sample_id_prefix}-{plate}-{multiplex_group}-{barcode_name}'

        # THIS IS BASED ON FORMAT BCL2FASTQ NEEDS
        records.append({'Sample_ID': sample_id,
                        'index': i7_barcode,
                        'index2': reverse_comp(i5_barcode) if i5_reverse_comp else i5_barcode,
                        'Sample_Project': row['tube_label'],
                        'Description': row['email']})

    miseq_sample_sheet, nova_sample_sheet = _make_final_samplesheet(records)
    return miseq_sample_sheet, nova_sample_sheet


def _make_final_samplesheet(records):
    # THIS IS BASED ON FORMAT BCL2FASTQ NEEDS
    sample_sheet = pd.DataFrame(records)
    sample_sheet['Sample_Name'] = ''
    sample_sheet['Sample_Well'] = ''
    sample_sheet['I7_Index_ID'] = ''
    sample_sheet['I5_Index_ID'] = ''
    sample_sheet['Sample_Plate'] = 'Plate'

    miseq_sample_sheet = sample_sheet[['Sample_ID', 'Sample_Name', 'Sample_Plate',
                                       'Sample_Well', 'I7_Index_ID', 'index',
                                       'I5_Index_ID', 'index2', 'Sample_Project',
                                       'Description']].copy()

    lanes = []
    for i in range(1, 5):
        lane_df = miseq_sample_sheet.copy()
        lane_df['Lane'] = i
        lanes.append(lane_df)
    nova_sample_sheet = pd.concat(lanes)
    nova_sample_sheet = nova_sample_sheet[['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                                           'Sample_Well', 'I7_Index_ID', 'index',
                                           'I5_Index_ID', 'index2', 'Sample_Project',
                                           'Description']].copy()
    return miseq_sample_sheet, nova_sample_sheet


def make_sample_sheet(plate_info_path: str, output_prefix: str, header_path=None):
    """
    Make two sample sheets for bcl2fastq based on the plate info file:
    one for miseq, one for novaseq.

    Parameters
    ----------
    plate_info_path
    output_prefix
    header_path
    """
    critical_info, plate_info = _read_plate_info(plate_info_path)

    for k in ['input_plate_size']:
        if critical_info[k] not in LIMITED_CHOICES[k]:
            raise ValueError(f'Invalid value in critical_info section for {k}: {critical_info[k]}')

    barcode_table = pd.read_csv(BARCODE_TABLE_PATH, sep='\t', index_col='set_384_plate_pos')

    for i5_reverse_comp in [True, False]:
        miseq_sample_sheet, nova_sample_sheet = _plate_384_random_index_384(
            plate_info, barcode_table, i5_reverse_comp=i5_reverse_comp)

        # write miseq sample sheet (only standard i5)
        if not i5_reverse_comp:
            with open(output_prefix + '.miseq.sample_sheet.csv', 'w') as output_f:
                if header_path is None:
                    output_f.write(SAMPLESHEET_DEFAULT_HEADER)
                else:
                    with open(header_path) as hf:
                        output_f.write(hf.read())
                output_f.write(miseq_sample_sheet.to_csv(index=None))

        # write novaseq sample sheet
        with open(output_prefix + f'.novaseq.{"V1.5" if i5_reverse_comp else "v1"}.sample_sheet.csv', 'w') as output_f:
            if header_path is None:
                output_f.write(SAMPLESHEET_DEFAULT_HEADER)
            else:
                with open(header_path) as hf:
                    output_f.write(hf.read())
            output_f.write(nova_sample_sheet.to_csv(index=None))
    return


def print_plate_info():
    """Print the V2 plate info template to stdout."""
    with open(PACKAGE_DIR / 'files/plate_info_template_v2.txt') as f:
        print(f.read())
