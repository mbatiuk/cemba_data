import pathlib

# Load defaults
CONFIG_DIR = pathlib.Path(__file__).parent / 'default_config'

MAPPING_MODE_CHOICES = ['mct','mct-multi', 'mc', 'm3c','m3c-multi', 'mc-multi']

def print_default_mapping_config(mode,
                                 genome_fasta,
                                 hisat3n_dna_ref,
                                 chrom_size_path,
                                 hisat3n_rna_ref=None,
                                 gtf=None,
                                 nome=False,
                                 **kwargs):
        mode = mode.lower()
        if mode.split('-')[0] not in MAPPING_MODE_CHOICES:
                raise ValueError(f'Unknown mode {mode}')

        # validate mct-only parameters
        if mode.split('-')[0] == 'mct':
                if hisat3n_rna_ref is None:
                        raise ValueError('hisat3n_rna_ref is required if mode is mct.')
                if gtf is None:
                        raise ValueError('gtf must be provided when mode is mct.')

        if mode.split('-')[0] == 'mc':
                if nome:
                        config_path = CONFIG_DIR / 'mapping_config_nome.ini'
                else:
                        config_path = CONFIG_DIR / 'mapping_config_mc.ini'
                with open(config_path) as f:
                        config_content = f.read()
        elif mode.split('-')[0] == 'mct':
                if nome:
                        config_path = CONFIG_DIR / 'mapping_config_mct-nome.ini'
                else:
                        config_path = CONFIG_DIR / 'mapping_config_mct.ini'
                with open(config_path) as f:
                        config_content = f.read()
                config_content = config_content.replace('CHANGE_THIS_TO_YOUR_HISAT3N_RNA_REFERENCE',
                                                        str(hisat3n_rna_ref))
                config_content = config_content.replace('CHANGE_THIS_TO_YOUR_GENE_ANNOTATION_GTF', str(gtf))
        elif mode.split('-')[0] == 'm3c':
                config_path = CONFIG_DIR / 'mapping_config_m3c.ini'
                with open(config_path) as f:
                        config_content = f.read()
        else:
                raise ValueError(f'Unknown mode {mode}')

        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_MODE', mode)
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_CHROM_SIZE_PATH', str(chrom_size_path))
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_HISAT3N_DNA_REFERENCE', str(hisat3n_dna_ref))
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_REFERENCE_FASTA', str(genome_fasta))
        print(config_content)
        for key in kwargs:
                if kwargs[key] is None:
                        continue
                print(f"{key} = {kwargs[key]}")
        return
