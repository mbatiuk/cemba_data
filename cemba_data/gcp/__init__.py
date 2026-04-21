from .yap_gcp import *
from ..mapping.pipelines.mhap import bam2mhap
import fire

def main():
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire({
		'run_demultiplex': run_demultiplex,
		'run_mapping': run_mapping,
		'start_from_cell_bam': start_from_cell_bam,
		'bam2mhap': bam2mhap,
	}, serialize=lambda x: print(x) if not x is None else print(""))


if __name__ == "__main__":
	main()