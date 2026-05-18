import configparser
import pathlib
import pandas as pd
from papermill import execute_notebook, PapermillExecutionError


MODULE_DIR = pathlib.Path(__file__).parent

# plate info
def _parse_cell_id(cell_id):
    plate, multiplex_group, random_index = cell_id.split('-')
    # 384 pos
    col384 = int(random_index[1:]) - 1
    row384 = ord(random_index[0]) - 65  # convert A-P to 0-23
    return pd.Series({
        'Plate': plate,
        'MultiplexGroup': multiplex_group,
        'RandomIndex': random_index,
        'Col384': col384,
        'Row384': row384
    })


def get_plate_info(cell_ids):
    try:
        plate_info = pd.DataFrame(
            [_parse_cell_id(cell_id) for cell_id in cell_ids], index=cell_ids
        )
    except Exception:
        print(
            'Errors occur during parsing the plate info, this happens '
            'when the input FASTQ file name is not following naming convention. '
            'The `yap summary` also cannot generate html report due to missing plate info. '
            'In this case, you need to add the plate info by yourself to make plate view plots. '
            'This information is not necessary for following analysis though.'
        )
        plate_info = pd.DataFrame([], index=cell_ids)
    return plate_info


# Final summary function. It will aggregate all mapping summaries
# It will also add Plate, MultiplexGroup, RandomIndex, Col384, and Row384 metadata
# And it will also write ALLC path file for generating MCDS

def final_summary(output_dir, notebook=None, mode=None, kernel_name='python3'):
    output_dir = pathlib.Path(output_dir).absolute()

    # Auto-detect mode from the mapping config saved in the output directory
    if mode is None:
        config_path = output_dir / 'mapping_config.ini'
        if config_path.exists():
            cfg = configparser.ConfigParser()
            cfg.read(config_path)
            mode = cfg.get('mode', 'mode', fallback=None)
        if mode is None:
            raise ValueError(
                'Could not determine mode automatically. '
                'Please specify --mode (e.g. mc, mct, or m3c).'
            )

    mode = mode.split('-')[0]

    # Ensure all UID dirs with a Snakefile also have MappingSummary (successful run)
    snakefile_list = list(output_dir.glob('*/Snakefile'))
    summary_paths = []
    missing_summary_dirs = []
    for path in snakefile_list:
        uid_dir = path.parent
        summary_path = uid_dir / 'MappingSummary.csv.gz'
        if summary_path.exists():
            summary_paths.append(summary_path)
        else:
            missing_summary_dirs.append(uid_dir)

    if len(missing_summary_dirs) != 0:
        print('These sub dirs are missing MappingSummary files:')
        for p in missing_summary_dirs:
            print(p)
        raise FileNotFoundError(
            f'All sub dirs must be successfully mapped before generating final summary.\n'
            f'MappingSummary.csv.gz is the final target of the snakefile.\n'
            f'Run the corresponding snakemake command again to retry mapping.\n'
            f'Commands can be found in output_dir/snakemake/*/snakemake_cmd.txt'
        )

    # aggregate UID-level mapping summaries
    total_mapping_summary = pd.concat(
        [pd.read_csv(path, index_col=0) for path in summary_paths]
    )
    total_mapping_summary_path = output_dir / 'stats/MappingSummary.csv.gz'

    # plate info
    _plate_info = get_plate_info(total_mapping_summary.index)
    total_mapping_summary = pd.concat([_plate_info, total_mapping_summary], axis=1)

    # save aggregated summary
    total_mapping_summary.to_csv(total_mapping_summary_path)

    # write ALLC path file for generating MCDS
    allc_paths = pd.Series({
        path.name.split('.')[0]: str(path)
        for path in output_dir.glob('*/allc/*tsv.gz')
    })
    allc_paths.to_csv(output_dir / 'stats/AllcPaths.tsv', sep='\t', header=False)

    if 'Plate' in total_mapping_summary.columns:  # only run notebook when plate info exists
        nb_path = output_dir / 'stats/MappingSummary.ipynb'
        try:
            if notebook is None:
                template_notebook = MODULE_DIR / f'mapping_summary_template/{mode}_template.ipynb'
            else:
                template_notebook = str(notebook)

            print(f'Using notebook template from {template_notebook}')
            print(f'Executing summary plotting notebook using kernel {kernel_name}')
            execute_notebook(
                input_path=str(template_notebook),
                output_path=str(nb_path),
                parameters=dict(output_dir=str(output_dir)),
                kernel_name=kernel_name,
            )
            print('Summary notebook successfully executed. Exporting HTML...')
            import subprocess
            subprocess.run(['jupyter', 'nbconvert', '--to', 'html', str(nb_path)])
            print(f'See the summary plots here: {str(nb_path)[:-5]}html')
            print(f'Or customise the summary plots here: {nb_path}')
        except PapermillExecutionError:
            print(f'Summary plotting notebook encountered an error. Check: {nb_path}')
    return


