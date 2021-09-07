import json
from yaml import safe_load

from definitions import CONF_PATH


def get_conf():
    """Read the config.yaml file for the project"""
    with open(CONF_PATH, 'r') as file:
        conf = safe_load(file)
    return conf


def get_product_dict(exp_dir):
    """
    Read the compound_alternative_mass_dict.json and return the dictionary linking SumFxx and A,B,C... product
    identifiers.
    :param exp_dir: directory or the experiment under consideration
    :return: dict of the form {"A_formula": "SumF1", ...}
    """
    with open(exp_dir / 'compound_alternative_mass_dict.json', 'r') as file:
        product_dict = json.load(file)
    return product_dict


def get_internal_standard_number(exp_dir):
    """
    Read the SumFxx number of the internal standard for an experiment
    :param exp_dir: directory or the experiment under consideration
    :return: str, SumFxx number for the internal standard
    """
    product_dict = get_product_dict(exp_dir)
    return product_dict['IS_formula']
