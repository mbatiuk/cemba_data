import collections
import functools
from typing import List


@functools.lru_cache(maxsize=10)
def parse_chrom_size(path, remove_chr_list=None):
    """
    Parse UCSC chrom size file.

    Support simple UCSC chrom size file, or .fai format (1st and 2nd columns same as chrom size file)
    """
    if remove_chr_list is None:
        remove_chr_list = []

    with open(path) as f:
        chrom_dict = collections.OrderedDict()
        for line in f:
            # *_ for other format like fadix file
            chrom, length, *_ = line.strip("\n").split("\t")
            if chrom in remove_chr_list:
                continue
            chrom_dict[chrom] = int(length)
    return chrom_dict



def genome_region_chunks(chrom_size_path: str, bin_length: int = 10000000, combine_small: bool = True) -> List[str]:
    """
    Split the whole genome into bins, where each bin is {bin_length} bp. Used for tabix region query.

    Parameters
    ----------
    chrom_size_path
        Path of UCSC genome size file
    bin_length
        length of each bin
    combine_small
        whether combine small regions into one record

    Returns
    -------
    list of records in tabix query format
    """
    chrom_size_dict = parse_chrom_size(chrom_size_path)

    cur_chrom_pos = 0
    records = []
    record_lengths = []
    for chrom, chrom_length in chrom_size_dict.items():
        while cur_chrom_pos + bin_length <= chrom_length:
            # tabix region is 1 based and inclusive
            records.append(f"{chrom}:{cur_chrom_pos}-{cur_chrom_pos + bin_length - 1}")
            cur_chrom_pos += bin_length
            record_lengths.append(bin_length)
        else:
            records.append(f"{chrom}:{cur_chrom_pos}-{chrom_length}")
            record_lengths.append(chrom_length - cur_chrom_pos)
            cur_chrom_pos = 0

    # merge small records (when bin larger then chrom length)
    final_records = []
    if combine_small:
        temp_records = []
        cum_length = 0
        for record, record_length in zip(records, record_lengths):
            temp_records.append(record)
            cum_length += record_length
            if cum_length >= bin_length:
                final_records.append(" ".join(temp_records))
                temp_records = []
                cum_length = 0
        if len(temp_records) != 0:
            final_records.append(" ".join(temp_records))
    else:
        for record in records:
            final_records.append(record)
    return final_records




