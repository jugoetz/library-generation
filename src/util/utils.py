import json
from yaml import safe_load
import pandas as pd

from src.definitions import CONF_PATH, PLATE_LIST_PATH


def get_conf():
    """Read the config.yaml file for the project"""
    with open(CONF_PATH, "r") as file:
        conf = safe_load(file)
    return conf


def get_product_dict(exp_dir):
    """
    Read the compound_alternative_mass_dict.json and return the dictionary linking SumFxx and A,B,C... product
    identifiers.

    :param exp_dir: directory or the experiment under consideration
    :return: dict of the form {"A_formula": "SumF1", ...}
    """
    with open(exp_dir / "compound_alternative_mass_dict.json", "r") as file:
        product_dict = json.load(file)
    return product_dict


def get_internal_standard_number(exp_dir):
    """
    Read the SumFxx number of the internal standard for an experiment

    :param exp_dir: Directory of the experiment under consideration
    :type exp_dir: pathlib.Path
    :return: SumFxx number for the internal standard
    :rtype: str
    """
    product_dict = get_product_dict(exp_dir)
    return product_dict["IS_formula"]


def get_lcms_file_name(lab_journal_number):
    """
    :param lab_journal_number: JGxxx, the lab journal number attached to the experiment
    :type lab_journal_number: str
    :return: filename of the lcms results file
    :rtype: str
    """
    df = pd.read_csv(PLATE_LIST_PATH)
    return df.loc[
        df["lab_journal_number"] == lab_journal_number, "results_file_name"
    ].values[0]
