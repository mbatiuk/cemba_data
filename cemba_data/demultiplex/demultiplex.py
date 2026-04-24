"""
Demultiplex pipeline
"""

import logging
import re
import pandas as pd

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


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

