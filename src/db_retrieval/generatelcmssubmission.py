"""
Generate the submission file for high-throughput MS of library compounds.

Inputs:
  - Plate layout (csv-files, usually 'plate_layout_plateX.csv'):
    One plate per file. Handles arbitrary number of plates/files. X may be any number.
  - EITHER
     the mol_prop_dict prepared from sdf_to_properties
    OR
     a SQLite DB table of the virtual library
  - Identity of compounds (txt-file):
    Space-delimited file with one compound per row.
    Columns contain short name (e.g. 'M1') used in the plate layout file
    and the long name (e.g. 2-Pyr003_MT or 2-Pyr003) used in the SDF files

Output:
  - Submission file (output.csv):
    One well per row.
    Columns are:
      1 well identifier (e.g. Py-A-3)
    + 8 molecular formulae for SynFerm products A-H
    + a variable number of deprotection product columns (e.g. from Boc removal)
    + 1 column for internal standard

TODO this is in dire need of refactoring
"""

import gzip
import json
import re
import os
import sqlite3 as sql
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd

from definitions import LIB_STATIC_DIR, PLATES_DIR, COMPOUND_MAPPING_PATH, CONF_PATH, DB_PATH
from labware.plates import Plate384, Plate96
from utils import get_conf

"""GLOBALS"""
PROP_SOURCE = 'db'  # ['db' / 'dict'] where to look up molecular formula or mass.
# dict only gives expected product props, db also has deprotection products
VERBOSE = True
ADD_IS = True
MASS_OR_FORMULA = 'formula'  # ['mass'/'formula'] can output mass as a number or can give chemical formula
PLATE_SIZE = 384
conf = get_conf()


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
    if conf['well_plate_size'] == 96:
        p = Plate96(150000, 5000)
    elif conf['well_plate_size'] == 384:
        p = Plate384(12000, 2500)
    else:
        raise ValueError(f'Invalid plate size: {conf["well_plate_size"]}')
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
    """reorder mol_prop_dict so that letters are in the outer level and long_names in the inner. Lets keep this out 
    for now and try to find a better way """
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
    df.loc[df['long'] != '', 'IS_mass'] = conf['is_mass']
    df.loc[df['long'] != '', 'IS_formula'] = conf['is_formula']
    return df


def write_csv(df, file, exp_dir):
    """
    Generate formatted csv output for Mobias
    """

    # TODO this may only take one plate (it does, but need to change how it is called)
    def splitwell(df):
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
    subset_df = df.loc[df['long'] != '', subset]  # remove columns with no long name set (those are emtpy wells e.g. A1)
    # define rename dict for column renaming.
    # should always map product A_x -> SumF1 ... H_x -> SumF8, IS -> SumF9
    rename_dict = {}
    counter = 1
    for column in subset_df.columns:
        if column.endswith(suffix):
            rename_dict[column] = f'SumF{counter}'
            counter += 1

    subset_df.rename(columns=rename_dict, inplace=True)
    with open(exp_dir / 'compound_alternative_mass_dict.json', 'w') as jsonfile:
        json.dump(rename_dict, jsonfile)  # save this

    subset_df.replace('n/a', '',
                      inplace=True)  # we have n/a strings where enumeration has not given a product. Replace them with empty string.
    # finally, resort the dataframe columns so that vial is in front and SumFx ascending
    cols = subset_df.columns.tolist()
    cols.remove('Vial')
    reordered_columns = ['Vial'] + sorted(cols, key=lambda x: int(x.strip('SumF')))
    subset_df = subset_df[reordered_columns]
    subset_df.to_csv(file, index=False)
    return


def main(exp_dir):
    if PROP_SOURCE == 'dict':
        """Import Mass and Formula from json"""
        with gzip.open(LIB_STATIC_DIR / 'static_mol_prop_dict.json.gz', 'rt', encoding='ascii') as zipfile:
            mol_prop_dict = json.load(zipfile)
    elif PROP_SOURCE == 'db':
        """Import Mass and Formula from DB"""
        mol_prop_dict = get_prop_from_db(DB_PATH)
    else:
        raise ValueError(f'Invalid PROP_SOURCE {PROP_SOURCE}')

    """Import identities"""
    starting_material_dict = import_sm(COMPOUND_MAPPING_PATH)

    """Import plates as lists"""
    plates_dict = {}
    for path, _, files in os.walk(exp_dir):
        for f in files:
            m = re.compile(conf['plate_regex']).match(f)
            if m:
                plate_dict = import_pl(Path(path, f))
                plates_dict[m.group(1)] = plate_dict

    if VERBOSE:
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

    if VERBOSE:
        for i, data in df.iterrows():
            print(f'Plate {data["plate"]}, Well {data["well"]}, Product: {data["long"]}')

    """
    Add molecular formulae and exact masses to dataframe
    mol_prop_dict looks a little different depending on the source, so for now, we will go with two ways to read it
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

        combinations_a = len(set([key for v in needed_dict.values() for key in v.keys() if key.startswith('A')]))
        combinations_d = len(set([key for v in needed_dict.values() for key in v.keys() if key.startswith('D')]))
        combinations = max(combinations_a, combinations_d)
        # now we add everything to dataframes
        dfs = []
        # here we generate all those de-PG variants of the expected products
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

    """
    from the individual dataframes for different -PG product sets, we form one big dataframe and drop the na columns
    """
    df = dfs[0]
    for i in dfs[1:]:
        df = pd.merge(df, i, how='left')

    """add internal standard if user wishes"""
    if ADD_IS is True:
        add_is(df)
    else:
        print("You chose not to add internal standard.\n")

    """save to file"""
    output_file = exp_dir / 'mobias_submission.csv'
    write_csv(df, output_file, exp_dir)
    print(f'Data was written to "{output_file}".')

    return


if __name__ == '__main__':
    exp_nrs = ['JG264',
               'JG265',
               'JG266',
               'JG267',
               'JG268',
               'JG269',
               ]
    for exp_nr in exp_nrs:
        print(f'Generating submission file for {exp_nr}...')
        main(PLATES_DIR / exp_nr)
    print('End of script. Exiting...')
