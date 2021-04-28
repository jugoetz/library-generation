"""
Generate the submission file for high-throughput MS of library compounds.

Inputs:
  - Plate layout (csv-files, usually 'plate_layout_plateX.csv'):
    One plate per file. Handles arbitrary number of plates/files. X may be any number.
  - Compressed product SDF files (sdf.gz-files):
    Must contain long name and mol object
    RECENTLY CHANGED TO
    either the mol_prop_dict prepared from sdf_to_properties
        OR
    an SQLite DB of the virtuallibrary

  - Identity of compounds (txt-file):
    Space-delimited file with one compound per row.
    Columns contain short name (e.g. 'M1') used in the plate layout file
    and the long name (e.g. 2-Pyr003_MT or 2-Pyr003) used in the SDF files



Output:
  - Submission file (output.csv):
    One well per row. Columns are the well identifier (e.g. 1-A,3)
    up to six molecular formulae for SynFerm product_generator A-E + internal standard
"""

import pandas as pd
import re
import os
from pathlib import Path
from labware.plates import Plate384, Plate96
import pickle as pkl
import json
import gzip
import numpy as np
import sqlite3 as sql
from copy import deepcopy
from config import *

"""GLOBALS"""
# debug
USE_PICKLED_DF = False  # this will skip most of the script (if it has been run before). Only for debugging
# output controls
PROP_SOURCE = 'db'  # ['db' / 'dict'] where to look up molecular formula or mass. dict only gives expected product props, db also has deprotection products


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


def import_pl(file):
    """
    Import the plate layout from x files 'plate_layout_plateX.csv' where X is the plate number.
    Create a dict of the form {'A1': ['I1', 'M1', 'T1']} that maps wells to shorthand names
    :return: dict
    """
    if PLATE_SIZE == 96:
        p = Plate96(150000, 5000)
    elif PLATE_SIZE == 384:
        p = Plate384(12000, 2500)
    else:
        raise ValueError(f'Invalid plate size: {PLATE_SIZE}')
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
                or row[col] == 'None' \
                or row[col] == '' \
                or row[col] == 'n/a':
            pass  # catch empty fields (this will usually happen)
        else:
            try:
                long.append(dictionary[row[col]])
            except KeyError:  # catch unknown shorthand names (shouldn't usually happen)
                print(f'WARNING: Building block {row[col]} not found')
                long.append('n/a')
    return ' + '.join(long)


def get_prop_from_db(dbpath):
    """
    Get mass and formula from sqlite database. If there are additional formulae/masses in a designated column with _alt
    suffix, get those, too.
    Return a dictionary of all compounds in the db
    Final dictionary structure:
    {'Al036 + Mon003 + TerTH011: {'A': ('C5H8O', 210.21321), 'B'...}}
    if additional masses are present in DB:
    {'Al036 + Mon003 + TerTH011: {'A': ('C5H8O', 210.21321), 'B'..., 'A_2': ('C3H6O',150.3525), 'A_3': ('C3H5O',140.1251),}}
    consequently these will always be true:
    len(prop_dict) = number of unique long_names in db
    len(prop_dict[X]) = 8 if no protecting group, else 8 + number of deprotection combinations
    :param dbpath: str or path-like
    :return: dict
    """
    con = sql.connect(dbpath)
    cur = con.cursor()
    prop_dict = {}
    for row in cur.execute(
            'SELECT id, long_name, type, molecular_formula_1, molecular_formula_alt, lcms_mass_1, lcms_mass_alt FROM virtuallibrary').fetchall():
        if not row[1] in prop_dict:  # if we have not encountered this long name before, make a dictionary as the entry
            prop_dict[row[1]] = {}
        prop_dict[row[1]][row[2]] = (row[3], row[5])  # these are the "standard" formula and mass
        if row[4] is not None:  # if alternate masses are present
            for i, (f, m) in enumerate(zip(row[4].split(','), row[6].split(','))):
                prop_dict[row[1]][f'{row[2]}_{i + 2}'] = (f, m)
    con.close()
    """reorder mol_prop_dict so that letters are in the outer level and long_names in the inner. Lets keep this out for now and try to find a better way"""
    # m_prop_dict = {}
    # for long, val in prop_dict.items():
    #     for letter, value in val.items():
    #         try:
    #             m_prop_dict[letter][long] = value
    #         except KeyError:
    #             m_prop_dict[letter] = {long: value}
    return prop_dict


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


def add_is(df):
    df['IS_mass'] = np.nan
    df['IS_formula'] = ''
    df.loc[df['long'] != '', 'IS_mass'] = IS_MASS
    df.loc[df['long'] != '', 'IS_formula'] = IS_FORMULA
    return df


def write_csv(df, file):
    """
    Generate formatted csv output for Mobias
    """

    # TODO this may only take one plate (it does, but need to change how it is called)
    def splitwell(df):
        plate = df['plate']
        well = df['well']
        new = f'Py-{str(well)[0]}-{str(well)[1:]}'
        return new

    if MASS_OR_FORMULA == 'mass':
        suffix = '_mass'
    elif MASS_OR_FORMULA == 'formula':
        suffix = '_formula'
    else:
        raise ValueError(f'Invalid option {MASS_OR_FORMULA}')
    df['Vial'] = df.loc[:, ['plate', 'well']].apply(splitwell, axis=1)
    # drop any row were all masses are np.nan
    df.dropna(axis=0, how='all', subset=[col for col in df.columns if col.endswith('mass') and col != 'IS_mass'],
              inplace=True)
    subset = ['Vial', ]
    for s in df.columns:
        if s.endswith(suffix):
            subset.append(s)
    subset_df = df.loc[df['long'] != '', subset]  # remove columns with no long name set
    # define rename dict for column renaming.
    # should always map product A_x -> SumF1 ... H_x -> SumF8, IS -> SumF9
    rename_dict = {}
    for column in subset_df.columns:
        if column.startswith('A'):
            rename_dict[column] = 'SumF1'
        elif column.startswith('B'):
            rename_dict[column] = 'SumF2'
        elif column.startswith('C'):
            rename_dict[column] = 'SumF3'
        elif column.startswith('D'):
            rename_dict[column] = 'SumF4'
        elif column.startswith('E'):
            rename_dict[column] = 'SumF5'
        elif column.startswith('F'):
            rename_dict[column] = 'SumF6'
        elif column.startswith('G'):
            rename_dict[column] = 'SumF7'
        elif column.startswith('H'):
            rename_dict[column] = 'SumF8'
        elif column.startswith('IS'):
            rename_dict[column] = 'SumF9'

    subset_df.rename(columns=rename_dict, inplace=True)
    # sometimes we might not have all columns (since the products don't exist). In that case add the missing ones and fill with empty values
    necessary_columns = [f'SumF{i + 1}' for i in range(9)]
    for col in necessary_columns:
        if col not in subset_df.columns:
            subset_df[col] = ''
    subset_df.replace('n/a', '',
                      inplace=True)  # we have n/a strings where enumeration has not given a product. Replace them with empty string.
    # finally, resort the dataframe columns
    subset_df = subset_df[['Vial', 'SumF1', 'SumF2', 'SumF3', 'SumF4', 'SumF5', 'SumF6', 'SumF7', 'SumF8', 'SumF9', ]]
    subset_df.to_csv(file, index=False)
    return


########################################################################################################################
#                                                        MAIN                                                          #
########################################################################################################################


if __name__ == '__main__':
    if USE_PICKLED_DF is False:
        if PROP_SOURCE == 'dict':
            """Import Mass and Formula from json"""
            with gzip.open(SDF_DIR / 'static_mol_prop_dict.json.gz', 'rt', encoding='ascii') as zipfile:
                mol_prop_dict = json.load(zipfile)
        elif PROP_SOURCE == 'db':
            """Import Mass and Formula from DB"""
            mol_prop_dict = get_prop_from_db(DB_PATH)
        else:
            raise ValueError(f'Invalid PROP_SOURCE {PROP_SOURCE}')

        """Import identities"""
        starting_material_dict = import_sm(COMPOUND_MAPPING)

        """Import plates as lists"""
        plates_dict = {}
        for path, _, files in os.walk(EXP_DIR):
            for f in files:
                m = PLATE_REGEX.match(f)
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

        """
        Add molecular formulae and exact masses to dataframe
        mol_prop_dict looks a little different dependign on the source so for now, we will go with two ways to read it
        """
        if PROP_SOURCE == 'dict':
            for letter, mol_props in sorted(mol_prop_dict.items()):
                df[f'{letter}_mass'] = df['long'].apply(get_prop, mol_props=mol_props, prop='mass')
                df[f'{letter}_formula'] = df['long'].apply(get_prop, mol_props=mol_props, prop='formula')
        elif PROP_SOURCE == 'db':
            # select the subportion of the prop-dict that we need (for this specific plate)
            needed_dict = {k: v for k, v in mol_prop_dict.items() if k in df['long'].values}
            # find how many different submission forms we will need
            # Here's how this works: We only need look at products A and E since they are candidates for the
            # maximum number of combinations. Everything else will have equal or less.
            # We go through the keys in the inner dict-level, only take unique ones starting with A, resp. D. The length
            # of this set corresponds to the maximum number of combinations and thus submission forms

            combinations_A = len(set([key for v in needed_dict.values() for key in v.keys() if key.startswith('A')]))
            combinations_D = len(set([key for v in needed_dict.values() for key in v.keys() if key.startswith('D')]))
            combinations = max(combinations_A, combinations_D)
            # now we add everything to dataframes
            dfs = []
            for i in range(combinations):
                index = i + 1
                dfs.append(deepcopy(df))
                df_this = dfs[i]
                for long_name, type_dict in needed_dict.items():
                    for letter, mol_props in type_dict.items():
                        if index == 1 and letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
                            # if the columns are not yet present, fill them with placeholder values
                            if f'{letter}_mass' not in df_this.columns:
                                df_this[f'{letter}_mass'] = np.nan
                            if f'{letter}_formula' not in df_this.columns:
                                df_this[f'{letter}_formula'] = ''
                            df_this.loc[df_this['long'] == long_name, f'{letter}_mass'] = mol_props[1]
                            df_this.loc[df_this['long'] == long_name, f'{letter}_formula'] = mol_props[0]
                        if letter.endswith(str(index)):
                            # if the columns are not yet present, fill them with placeholder values
                            if f'{letter}_mass' not in df_this.columns:
                                df_this[f'{letter}_mass'] = np.nan
                            if f'{letter}_formula' not in df_this.columns:
                                df_this[f'{letter}_formula'] = ''
                            df_this.loc[df_this['long'] == long_name, f'{letter}_mass'] = mol_props[1]
                            df_this.loc[df_this['long'] == long_name, f'{letter}_formula'] = mol_props[0]

        """add internal standard if user wishes"""
        if ADD_IS == "y":
            if PROP_SOURCE == 'dict':
                add_is(df)
            elif PROP_SOURCE == 'db':
                for df in dfs:
                    add_is(df)
        else:
            print("You chose not to add internal standard.\n")

        """Output to pickle"""
        if PROP_SOURCE == 'dict':
            with open(EXP_DIR / 'dataframe.pkl', 'wb') as file:
                pkl.dump(df, file)
        elif PROP_SOURCE == 'db':
            with open(EXP_DIR / 'dataframes.pkl', 'wb') as file:
                pkl.dump(dfs, file)
    if PROP_SOURCE == 'dict':
        with open(EXP_DIR / 'dataframe.pkl', 'rb') as file:
            df = pkl.load(file)
    elif PROP_SOURCE == 'db':
        with open(EXP_DIR / 'dataframes.pkl', 'rb') as file:
            dfs = pkl.load(file)
    if PROP_SOURCE == 'dict':
        output_file = EXP_DIR / 'mobias_submission.csv'
        write_csv(df, output_file)
        print(f'Data was written to "{output_file}".')
    elif PROP_SOURCE == 'db':
        for i, df in enumerate(dfs):
            output_file = EXP_DIR / f'mobias_submission_{i + 1}.csv'
            write_csv(df, output_file)
            print(f'Data was written to "{output_file}".')
    print('End of script. Exiting...')
