from yaml import safe_load

from definitions import CONF_PATH


def get_conf():
    """Read the config.yaml file for the project"""
    with open(CONF_PATH, 'r') as file:
        conf = safe_load(file)
    return conf
