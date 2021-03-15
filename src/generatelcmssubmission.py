"""
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
import re
import os
from pathlib import Path
import gzip
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcMolFormula
from labware.plates import Plate384
import pickle as pkl
import numpy as np

"""GLOBALS"""
DATA_DIR = Path('..', 'data').resolve()
OUTPUT_DIR = DATA_DIR / 'outputs'
INPUT_DIR = DATA_DIR / 'inputs'
EXP_DIR = OUTPUT_DIR / 'target_plates' / 'exp1'
verbose = True
USE_PICKLED_MOLPROPDICT = True
USE_PICKLED_DF = True
ADD_IS = 'y'


def import_sm(file):
    """
    Import identity of starting materials into a dictionary.
    Return dictionary that maps shorthand names (e.g. 'T1') to long names (e.g. 'TerTH001')
    :param file: str or Path-like
    :return: dict
    """
    dictionary = {}
    with open(file, 'r') as f:
        # iterate txt-file
        for line in f:
            # split columns of space-delimited file
            columns = line.strip('\n').split(sep=' ', maxsplit=1)
            # remove any "_XX" suffixes from long names
            dictionary[columns[0]] = columns[1].split("_")[0]
    return dictionary


def import_mol(file):
    """
    Import product mol from gzipped SDF. Return a generator that yields rdkit mol objects
    :param file: str or Path-like
    :return: generator
    """
    if str(file).endswith('.sdf.gz'):
        supply = ForwardSDMolSupplier(gzip.open(file))
    else:
        raise ValueError(f'Invalid file "{file}". User must supply .sdf.gz file.')
    return supply


def import_pl(file):
    """
    Import the plate layout from x files 'plate_layout_plateX.csv' where X is the plate number.
    Create a dict of the form {'A1': ['I1', 'M1', 'T1']} that maps wells to shorthand names
    :return: dict
    """
    p = Plate384(12000, 2500)
    p.from_csv(file)
    p_dict = p.to_dict()
    return p_dict


def get_long_name(row, dictionary):
    """
    For use with pandas.DataFrame.apply()
    Convert shorthand names located in different df columns (I, M, T) to one longhand name.
    """
    long = []
    for col in ['I', 'M', 'T']:
        if row[col] is None \
                or row[col] == 'None'\
                or row[col] == ''\
                or row[col] == 'n/a':
            pass  # catch empty fields (this will usually happen)
        else:
            try:
                long.append(dictionary[row[col]])
            except KeyError:  # catch unknown shorthand names (shouldn't usually happen)
                print(f'WARNING: Building block {row[col]} not found')
                long.append('n/a')
    return ' + '.join(long)


def get_prop(long_name, mol_props, prop):
    """
    Get mass and formula of a molecule identified by long_name from the mol_props dictionary.
    :param long_name: str
    :param mol_props: dict
    :param prop: str, either of "mass" or "formula"
    """
    if prop == 'mass':
        i = 0
    elif prop == 'formula':
        i = 1
    else:
        raise ValueError('')
    if long_name == '':
        return None
    else:
        try:
            mass = mol_props[long_name][i]
        except KeyError:
            mass = 'n/a'
    return mass


def write_csv(df, file):
    """
    Generate formatted csv output for Mobias
    """
    def splitwell(df):
        plate = df['plate']
        well = df['well']
        new = str(plate) + '-' + str(well)[0] + ',' + str(well)[1:]
        return new
    df['vial'] = df.loc[:, ['plate', 'well']].apply(splitwell, axis=1)
    subset = ['vial',]
    for s in df.columns:
        if s.endswith('_formula'):
            subset.append(s)
    subset_df = df.loc[df['long'] != '', subset]
    rename_dict = {}
    i = 1
    for column in subset_df.columns:
        if column.endswith('_formula'):
            rename_dict[column] = f'SumF{i}'
            i += 1
    subset_df.rename(columns=rename_dict, inplace=True)
    subset_df.to_csv(file, index=False)
    return


########################################################################################################################
#                                                        MAIN                                                          #
########################################################################################################################


if __name__ == '__main__':
    if USE_PICKLED_DF is False:
        if USE_PICKLED_MOLPROPDICT is False:
            """Generate mol suppliers to import from SDF"""
            mol_prop_dict = {}
            for path, _, files in os.walk(DATA_DIR / 'library_static'):
                for f in files:
                    if f.startswith('product_') and f.endswith('.sdf.gz'):
                        letter = f.split('_')[1].split('.')[0]
                        supplier = import_mol(Path(path, f))
                        mol_props = {}
                        for mol in supplier:
                            if mol is not None:
                                mol_props[mol.GetProp('_Name')] = [CalcExactMolWt(mol), CalcMolFormula(mol)]
                        mol_prop_dict[letter] = mol_props
            with open('debug_mol_prp_dict.pkl', 'wb') as file:
                pkl.dump(mol_prop_dict, file)

        with open('debug_mol_prp_dict.pkl', 'rb') as file:
            mol_prop_dict = pkl.load(file)

        """Import identities (this can stay as is)"""
        starting_material_dict = import_sm(OUTPUT_DIR / 'compound_mapping.txt')

        """Import plates as lists"""
        plates_dict = {}
        plate_regex = re.compile('plate_layout_plate([0-9]+).csv')
        for path, _, files in os.walk(EXP_DIR):
            for f in files:
                m = plate_regex.match(f)
                if m:
                    plate_dict = import_pl(Path(path, f))
                    plates_dict[m.group(1)] = plate_dict

        if verbose:
            # print the input values for double checking
            print("########## INPUT VALUES ###########\n")
            for k, v in mol_prop_dict.items():
                print(f"Products {k}:\n{v}\n")
            print(f"Used starting materials: \n{starting_material_dict}\n\n")
            for k, v in plates_dict.items():
                print(f"Plate layout {k}:\n{v}\n")

            # start printing terminal output
            print("########## OUTPUT ###########\n")

        """Put everything in df for ease of use"""
        dfs = []
        for plate, plate_content in plates_dict.items():
            df = pd.DataFrame.from_dict(plate_content, orient='index', columns=['I', 'M', 'T'])
            df.reset_index(drop=False, inplace=True)
            df.rename({'index': 'well'}, axis=1, inplace=True)
            df.insert(loc=0, column='plate', value=int(plate))
            dfs.append(df)
        df = pd.concat(dfs)
        df.sort_values(by=['plate'], inplace=True, kind='mergesort')  # mergesort
        df.reset_index(drop=True, inplace=True)

        """Translate product names from shorthand to long names (e.g. 'Al002 + Mon001 + TerTH010')"""
        df['long'] = df.loc[:, ['I', 'M', 'T']].apply(get_long_name, axis=1, dictionary=starting_material_dict)

        if verbose:
            for i, data in df.iterrows():
                print(f'Plate {data["plate"]}, Well {data["well"]}, Product: {data["long"]}')

        """Add molecular formulae and exact masses to dataframe"""
        for letter, mol_props in sorted(mol_prop_dict.items()):
            df[f'{letter}_mass'] = df['long'].apply(get_prop, mol_props=mol_props, prop='mass')
            df[f'{letter}_formula'] = df['long'].apply(get_prop, mol_props=mol_props, prop='formula')

        """add internal standard if user wishes"""
        internal_standard_y_n = ADD_IS
        internalstandard_formula = "C20H21O4Cl"  # molecular formula of Fenofibrat
        internalstandard_mass = 360.1128  # exact mass of Fenofibrat
        if internal_standard_y_n == "y":
            df['IS_mass'] = np.nan
            df['IS_formula'] = ''
            df.loc[df['long'] != '', 'IS_mass'] = internalstandard_mass
            df.loc[df['long'] != '', 'IS_formula'] = internalstandard_formula
        else:
            print("You chose not to add internal standard.\n")

        """Output to pickle"""
        with open(EXP_DIR / 'dataframe.pkl', 'wb') as file:
            pkl.dump(df, file)

    with open(EXP_DIR / 'dataframe.pkl', 'rb') as file:
        df = pkl.load(file)
    output_file = EXP_DIR / 'mobias_submission.csv'
    write_csv(df, output_file)
    print(f'Data was written to "{output_file}".')
    print('End of script. Exiting...')
