import pathlib
import cemba_data
import pandas as pd
from papermill import execute_notebook, PapermillExecutionError


PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


# plate info
def _parse_cell_id_v1(cell_id):
    plate1, plate2, pcr_index, random_index = cell_id.split('-')
    if random_index.upper() in {'AD001', 'AD002', 'AD004', 'AD006'}:
        plate = plate1
    else:
        plate = plate2
    # 96 pos
    col96 = int(pcr_index[1:]) - 1
    row96 = ord(pcr_index[0]) - 65  # convert A-H to 0-8
    # 384 pos
    ad_index_384_dict = {
        'AD001': (0, 0),
        'AD002': (0, 1),
        'AD004': (1, 0),
        'AD006': (1, 1),
        'AD007': (0, 0),
        'AD008': (0, 1),
        'AD010': (1, 0),
        'AD012': (1, 1)
    }
    col384 = 2 * col96 + ad_index_384_dict[random_index][0]
    row384 = 2 * row96 + ad_index_384_dict[random_index][1]
    return pd.Series({
        'Plate': plate,
        'PCRIndex': pcr_index,
        'RandomIndex': random_index,
        'Col384': col384,
        'Row384': row384
    })


def _parse_cell_id_v2(cell_id):
    plate, multiplex_group, pcr_index, random_index = cell_id.split('-')
    # 384 pos
    col384 = int(random_index[1:]) - 1
    row384 = ord(random_index[0]) - 65  # convert A-P to 0-23
    return pd.Series({
        'Plate': plate,
        'PCRIndex': pcr_index,
        'MultiplexGroup': multiplex_group,
        'RandomIndex': random_index,
        'Col384': col384,
        'Row384': row384
    })


def get_plate_info(cell_ids, barcode_version='V2'):
    func = _parse_cell_id_v1 if barcode_version == 'V1' else _parse_cell_id_v2
    try:
        plate_info = pd.DataFrame(
            [func(cell_id) for cell_id in cell_ids], index=cell_ids
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
# It will also add Plate, PCRIndex, MultiplexGroup, RandomIndex, Col384, and Row384 metadata
# And it will also write ALLC path file for generating MCDS
 
def final_summary(output_dir, notebook=None,
                  mode='m3c', barcode_version='V2', kernel_name='python3'):
    output_dir = pathlib.Path(output_dir).absolute()
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
    _plate_info = get_plate_info(total_mapping_summary.index, barcode_version=barcode_version)
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
                template_notebook = PACKAGE_DIR / f'files/mapping_summary_template/{mode}_template.ipynb'
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


