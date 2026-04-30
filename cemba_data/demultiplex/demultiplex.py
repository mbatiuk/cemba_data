"""
Demultiplex pipeline
"""

import re
import pandas as pd
import os


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


def _parse_index_fasta(fasta_path):
	records = {}
	with open(fasta_path) as f:
		key_line = True
		for line in f:
			if key_line:
				key = line.lstrip('>').rstrip('\n')
				key_line = False
			else:
				value = line.lstrip('^').rstrip('\n')
				records[key] = value
				key_line = True
	return records


def _read_cutadapt_result(stat_path):
	"""
	Parser of cutadapt output
	"""
	with open(stat_path) as f:
		p = re.compile(
			r"Sequence: .+; Type: .+; Length: \d+; Trimmed: \d+ times")
		series = []
		total_pairs = -1
		for line in f:
			if line.startswith('Total read pairs processed'):
				total_pairs = line.split(' ')[-1]
				total_pairs = int(''.join(re.compile(r'\d').findall(total_pairs)))
			m = p.search(line)
			if m is not None:
				result_dict = {}
				for i in m.group().split('; '):
					k, v = i.split(': ')
					result_dict[k] = v
				result_series = pd.Series(result_dict)
				series.append(result_series)
		total_df = pd.DataFrame(series)
		total_df['Trimmed'] = total_df['Trimmed'].apply(
			lambda c: c.split(' ')[0]).astype(int)
		total_df['TotalPair'] = total_pairs
		total_df['Ratio'] = total_df['Trimmed'] / total_pairs
	return total_df

