import configparser

def get_configuration(config_path):
    """
    Read .ini config file from given path
    """
    if isinstance(config_path, configparser.ConfigParser):
        return config_path
    ref_path_config = configparser.ConfigParser()
    ref_path_config.read(config_path)

    total_config = {}
    for name, section in ref_path_config.items():
        for k, v in section.items():
            total_config[k] = v
    return total_config


MAPPING_MODE_CHOICES = ['mct','mct-multi', 'mc', 'm3c','m3c-multi', 'mc-multi']
