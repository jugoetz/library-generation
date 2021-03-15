#!/usr/bin/env python

"""
DEPRECATION NOTICE:
This script has been superseeded and should not be used any more.



Generate the submission file for high-throughput MS of library compounds.

Inputs:
  - Plate layout (csv-files, 'plate_layout_plateX.csv'):
    One plate per file. Handles arbitrary number of plates/files. X must start at 1 and increase.
  - Up to five product files (csv-files, 'A.csv'..'E.csv'):
    1 product per row with ID and molecular formula as columns.
  - Identity of compounds (txt-file):
    Space-delimited file with one compound per row.
    Columns contain identifier (e.g. 'M1') used in the plate layout file
    and the long name (e.g. 2-Pyr003_MT or 2-Pyr003) used in the product files
  - Plate number (interactive):
    Usually int, but str is allowed

Output:
  - Submission file (output.csv):
    One well per row. Columns are the well identifier (e.g. 1-A,3)
    up to six molecular formulae for SynFerm product_generator A-E + internal standard

WARNING:
    The function remove_one_proton (which can be called optionally from prompt4) holds the hardcoded assumption
    that Schrödinger fails to neutralize product C which thus has one proton in excess.
"""


import pandas as pd
import csv
import sys
import os
from molmass import Formula
from pathlib import Path


def import_sm(file):
    """
    Import identity of starting materials into a dictionary.
    Return dictionary that maps shorthand names (e.g. 'T1') to long names (e.g. 'TerTH001')
    :return: dict
    """
    dictionary = {}
    with open(file) as f:
        # iterate txt-file
        for line in f:
            # split columns of space-delimited file
            columns = line.strip('\n').split(sep=' ', maxsplit=1)
            # remove any "_XX" suffixes from long names
            dictionary[columns[0]] = columns[1].split("_")[0]
    return dictionary


def import_mf(file):
    """
    Import product molecular formulae
    :param file: Path to csv file with 2 columns: product decriptor and molecular formula
    :return: dict
    """
    try:
        # import molecular formulae to pandas dataframe, then convert to dictionary
        mf_dict = pd.read_csv(file, index_col="row ID").to_dict()
        mf_dict = mf_dict["Molecular formula"]  # remove the outer level of the dictionary
        # strip the person identifiers and [1/1] from the titles, so that only say "2-Pyr002" remains
        new_dict = {}
        for dic_key in mf_dict.keys():
            old_key = dic_key.split(sep="+")
            new_key = []
            for element in old_key:
                new_key.append(element.split(sep="_")[0].strip().split(sep=" ")[0])
            new_key = " + ".join(new_key)
            new_dict[new_key] = mf_dict[dic_key]
    except FileNotFoundError:
        # triggered if csv file does not exist
        new_dict = {"undefined": ""}
    return new_dict


def import_pl():
    """
    Import the plate layout from x files 'plate_layout_plateX.csv' where X is the plate number.
    Create a dict of the form {"1-A,1": "I1, M1, T1"} that maps wells to shorthand names
    :return: dict
    """
    # let user specify number of plates
    number_of_plates = input("Please enter number of plates: \n")
    # import all plate layouts to single dictionary
    layout = {}
    # iterate plates
    for plate in range(1, int(number_of_plates)+1):
        filename = "plate_layout_plate{}.csv".format(plate)
        if verbose:
            print("Trying to open {}...\n".format(filename))
        rows = []
        columns = []
        wells = []
        i = True
        with open(filename, newline='') as csvfile:
            csvreader = csv.reader(csvfile, dialect="excel")
            for line in csvreader:
                if i:  # True for the first line
                    columns = line
                    i = False
                rows.append(line[0])  # Put the first element of every line (row letter) into list
                wells.append(line)  # put everything from csv file into list
        # prepare the dictionary
        for i in range(1, len(rows)):  # iterate excluding the first row (column numbers)
            for j in range(1, len(columns)):  # iterate excluding the first column (row letters)
                # construct key as e.g. "1-A,1" and assign value e.g. "I1, M1, T1"
                layout["".join([str(plate), "-", rows[i], ",", columns[j]])] = wells[i][j]
    return layout


def translate_names(building_block_string, SM_dictionary):
    """
    Generate the long name of reagents from shorthand "M1, I1, T1" string.
    Called separately for every well.
    :return: list
    """
    # Set n/a for unspecified wells (blanks)
    if building_block_string == "":
        return "n/a"
    # remove any whitespaces (to prevent errors) and split on commas
    building_blocks = building_block_string.replace(" ", "").split(sep=",")
    # get in alphabetic order
    building_blocks = sorted(building_blocks)
    # iterate the bueilding blocks
    for i in range(0, len(building_blocks)):
        # substitute shorthand by long names
        building_blocks[i] = SM_dictionary[building_blocks[i]]
    # return a list (length == # of BBs per well) with names of I, M, T in alphabetic order (of the SHORThand names)
    return building_blocks


def generate_product_names(building_blocks):
    """
    Concatenate the compound list to product long names
    :param building_blocks:
    :return:
    """
    if building_blocks == "n/a":
        return "n/a", "n/a", "n/a"
    product_name_abc = " + ".join(building_blocks)
    try:
        product_name_d = " + ".join([building_blocks[0], building_blocks[2]])
    except IndexError:
        product_name_d = product_name_abc
        product_name_abc = "undefined"
    try:
        product_name_e = building_blocks[2]
    except IndexError:
        product_name_e = "undefined"
    return product_name_abc, product_name_d, product_name_e


def remove_one_proton(dictionary, position):
    """
    Take a dictionary of the form {"key":[<molecular formula 1>, molecular formula 2>]...}.
    Remove one proton from the list element specified by position for every dict entry.
    :param dictionary:
    :param position:
    :return:
    """
    import re
    count = 0  # counts proton substractions
    for dic_key in dictionary.keys():
        str_list_to_modify = dictionary[dic_key][position].split(sep=" ")
        str_to_modify = str_list_to_modify[1]  # THIS HOLDS THE HARDCODED ASSUMPTION THAT H COMES SECOND IN THE MOLECULAR FORMULA!
        match = re.match(r"([a-z]+)([0-9]+)", str_to_modify, re.I)
        if match:
            items = match.groups()
            protons = "".join([items[0], str(int(items[1])-1)])
            str_list_to_modify[1] = protons
            dictionary[dic_key][position] = " ".join(str_list_to_modify)
            count += 1
        else:
            print("ERROR in function remove_one_proton: Failed to find proton count.")
            exit(1)
    if verbose:
        print("Removed a total of {} protons from {} product_generator.\n".format(count, len(dictionary.keys())))
        print("Updated dictionary after proton removal:\n{}\n".format(dictionary))
    return dictionary


def add_internal_standard(dictionary, molform):
    """
    Take a dictionary of the form {"key":[<molecular formula 1>, molecular formula 2>]...}.
    Append an entry with molform to every dictionary key
    :param dictionary:
    :param molform:
    :return:
    """
    for dic_key in dictionary.keys():
        dictionary[dic_key].append(molform)
    if verbose:
        print("Added internal standard.\n New molecular formulae table:\n{}\n".format(dictionary))
    return dictionary


def convert_mf_to_mass(dictionary):
    """
    Take a dictionary of the form {"key":[<molecular formula 1>, molecular formula 2>]...}
    Convert all molecular formulae to exact masses
    Return a dictionary of the form {"key":[<exact mass 1>, exact mass 2>]...}
    :param dictionary:
    :return:
    """
    for dic_key in dictionary.keys():
        for i in range(0, len(dictionary[dic_key])):
            if dictionary[dic_key][i] != '':  # omit empty entries
                # create a Formula object from string
                formula_object = Formula(dictionary[dic_key][i])
                # calculate monoisotopic mass
                dictionary[dic_key][i] = formula_object.isotope.mass
            else:
                continue  # do nothing if entry is empty
    return dictionary


def get_experiment_number():
    """
    Obtain the experiment number directly from the name of the enclosing folder or from interactive prompt
    :return:
    """
    cwd = os.getcwd()
    parent_folder = cwd.split("/")[-1]
    # ask user to confirm default or choose alternative name
    alternative_name = input("Experiment number was detected as '{}' from directory name. "
                             "Press ENTER to confirm or write a different experiment number...\n".format(parent_folder))
    if alternative_name != "":
        parent_folder = alternative_name
    experiment_number = parent_folder + "-"
    return experiment_number


def write_csv(dictionary, file, experiment_number):
    """
    Generate formatted csv output
    :param dictionary:
    :param file:
    :param experiment_number:
    :return:
    """
    header = ["Sample ID", "Vial"]
    max_length = 0
    for dic_key in dictionary.keys():
        if len(dictionary[dic_key]) > max_length:
            max_length = len(dictionary[dic_key])
    for i in range(1, max_length+1):
        header.append("SumF"+str(i))
    with open(file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile, dialect="excel")
        writer.writerow(header)
        i = 1
        for dic_key in dictionary.keys():
            writer.writerow([experiment_number+"{0:0=4d}".format(i)] + [dic_key] + dictionary[dic_key])
            i += 1
    return


########################################################################################################################
#                                                        MAIN                                                          #
########################################################################################################################
v = False
DATA_DIR = Path('data').resolve()
INPUT_DIR = DATA_DIR / 'inputs'

if __name__ == '__main__':
    if '-v' in sys.argv or v is True:
        verbose = True
    else:
        verbose = False

    # import data
    molecularFormulae_A = import_mf(INPUT_DIR / 'A.csv')
    molecularFormulae_B = import_mf(INPUT_DIR / 'B.csv')
    molecularFormulae_C = import_mf(INPUT_DIR / 'C.csv')
    molecularFormulae_D = import_mf(INPUT_DIR / 'D.csv')
    molecularFormulae_E = import_mf(INPUT_DIR / 'E.csv')
    startingMaterial_dict = import_sm(DATA_DIR / 'outputs' / 'compound_mapping.txt')
    plateLayout_dict = import_pl()

    if verbose:
        # print the input values for double checking
        print("########## INPUT VALUES ###########\n")
        print("Products A:\n{}\n".format(molecularFormulae_A))
        print("Products B:\n{}\n".format(molecularFormulae_B))
        print("Products C:\n{}\n".format(molecularFormulae_C))
        print("Products D:\n{}\n".format(molecularFormulae_D))
        print("Products E:\n{}\n".format(molecularFormulae_E))
        print("Used starting materials: \n{}\n\n".format(startingMaterial_dict))
        print("Plate layout:\n{}\n\n".format(plateLayout_dict))

        # start printing terminal output
        print("########## OUTPUT ###########\n")

    # dictionary that will hold the full compound descriptors (list of 5) for every well (keys=wells)
    translatedplateLayout = {}
    # for every well, generate the product names
    for key in plateLayout_dict.keys():
        starting_material = translate_names(plateLayout_dict[key], startingMaterial_dict)
        product_ABC, product_D, product_E = generate_product_names(starting_material)
        translatedplateLayout[key] = [product_ABC, product_ABC, product_ABC, product_D, product_E]
        if verbose:
            print("Well {}".format(key))
            print("Product ABC:\n{}".format(product_ABC))
            print("Product D:\n{}".format(product_D))
            print("Product E:\n{}\n".format(product_E))

    # dictionary that will hold the molecular formulae (list of 5) for every well (keys=wells)
    molecularFormulae_dict = {}
    # for every well, lookup molecular formula corresponding to product names
    for key in translatedplateLayout.keys():
        if "n/a" not in translatedplateLayout[key]:  # don't process well if n/a
            try:
                molecularFormulae_dict[key] = [molecularFormulae_A[translatedplateLayout[key][0]]]
            except KeyError:
                molecularFormulae_dict[key] = [""]
            try:
                molecularFormulae_dict[key].append(molecularFormulae_B[translatedplateLayout[key][1]])
            except KeyError:
                molecularFormulae_dict[key].append("")
            try:
                molecularFormulae_dict[key].append(molecularFormulae_C[translatedplateLayout[key][2]])
            except KeyError:
                molecularFormulae_dict[key].append("")
            try:
                molecularFormulae_dict[key].append(molecularFormulae_D[translatedplateLayout[key][3]])
            except KeyError:
                molecularFormulae_dict[key].append("")
            try:
                molecularFormulae_dict[key].append(molecularFormulae_E[translatedplateLayout[key][4]])
            except KeyError:
                molecularFormulae_dict[key].append("")
        else:
            if verbose:
                print("Well {} is empty. If well {} is non-empty in the input, an error has occurred.".format(key, key))
    if verbose:
        print("\nPrepared dictionary with molecular formulae by well:\n{}\n".format(molecularFormulae_dict))

    # remove one proton from product C if user wishes
    remove_proton_y_n = input("Should the mass of product C be reduced by one hydrogen? "
                              "(to correct the frequent Schrödinger error) [enter y for yes]...\n")
    if remove_proton_y_n == "y":
        molecularFormulae_dict = remove_one_proton(molecularFormulae_dict, 2)  # the 2 equals product C
    else:
        print("You chose to skip proton removal.\n")

    # add internal standard if user wishes
    internal_standard_y_n = input("Was Fenofibrat added as internal standard? [enter y for yes]...\n")
    internalStandardMF = "C20 H21 O4 Cl"  # molecular formula of Fenofibrat
    if internal_standard_y_n == "y":
        molecularFormulae_dict = add_internal_standard(molecularFormulae_dict, internalStandardMF)
    else:
        print("You chose not to add internal standard.\n")

    # convert molecular formula to exact mass if user wishes
    convert_to_mass = input("Convert molecular formulae to mass? [enter y for yes]...\n")
    if convert_to_mass == "y":
        molecularFormulae_dict = convert_mf_to_mass(molecularFormulae_dict)
        print("Molecular formulae were converted to exact masses.")
        if verbose:
            print("Updated dictionary:\n{}\n".format(molecularFormulae_dict))

    # generate output
    output_file = DATA_DIR / 'outputs' / 'mobias_submission.csv'
    ex_number = get_experiment_number()
    write_csv(molecularFormulae_dict, output_file, ex_number)
    print("Data was written to '{}'.".format(output_file))
    print("End of script. Exiting...")
