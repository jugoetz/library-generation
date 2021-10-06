"""
Populate the 'experiments' table of the database.
Location data is taken from plate_layout files in either experiment directories (e.g. exp2, standard usage) or
in individual directories with lab journal numbers (e.g. JG228, non-standard - add things not within 50k project scope).
Compound data is retrieved from 'virtuallibrary' table.

TODO
    This could have additional functionality. Currently I perform some steps of the experiment setup manually:
    - Make new directories for the six plates inside plates directory
    - Move expX/plate_layout_y* to respective new directory
    - Add the experiments at the end of plates_list.csv
    In addition, this should be refactored
"""

import re
import os
import sqlite3
from datetime import datetime
from pathlib import Path

import pandas as pd

from definitions import PLATES_DIR, DB_PATH
from db_retrieval.generatelcmssubmission import import_pl
from utils import get_conf

conf = get_conf()

# control variables
experiment_directory = 'exp9'
experiment_number = 9  # can be int or None if not a canonical (50 k) plate
synthesis_date = datetime(2021, 10, 5).timestamp()
labj_nr_dict = {
    '1': 'JG270',
    '2': 'JG271',
    '3': 'JG272',
    '4': 'JG273',
    '5': 'JG274',
    '6': 'JG275',
}


def get_plates_for_experiment(exp_dir):
    """
    Import the plate layouts for an experiment into a dictionary of plates.
    """
    plates_dict = {}
    for path, _, files in os.walk(PLATES_DIR / f'{exp_dir}'):
        for f in files:
            m = re.compile(conf['plate_regex']).match(f)
            if m:
                print(f'Retrieving plate layout {f}...')
                plate = import_pl(Path(path, f), return_type='plate')
                plates_dict[m.group(1)] = plate
    return plates_dict


def bulk_data_retrieval_from_virtuallibrary(con, initiators: list, monomers: list, terminators: list):
    """
    From the virtuallibrary table, retrieve all combinations of the all initiators,
    monomers and terminators given to the function.
    """
    cur = con.cursor()
    results = cur.execute(f'''SELECT id, initiator, monomer, terminator, initiator_long, monomer_long, terminator_long,
     long_name, type, SMILES FROM main.virtuallibrary 
     WHERE initiator IN ({(", ".join("?" for _ in initiators))}) 
     AND monomer IN ({(", ".join("?" for _ in monomers))}) 
     AND terminator IN ({(", ".join("?" for _ in terminators))});''',
                          initiators + monomers + terminators).fetchall()
    return results


def bulk_data_insertion_to_experiments(con, df):
    """
    Insert rows in the experiments table.
    Note that (as per the DB table layout), only one vl_id can be given for a row. We use vl_id_A.

    :param df: pandas.DataFrame with the columns ['vl_id_A', 'exp_nr', 'plate', 'well', 'lab_journal_nr', 'synthesis_date_unixepoch',
           'initiator', 'monomer', 'terminator', 'initiator_long', 'monomer_long', 'terminator_long', 'long_name', 'SMILES_A',
           'SMILES_B', 'SMILES_C', 'SMILES_D', 'SMILES_E', 'SMILES_F', 'SMILES_G', 'SMILES_H']
    :return: None
    """
    cur = con.cursor()
    cur.executemany(
        'INSERT INTO main.experiments (vl_id, exp_nr, plate_nr, well, lab_journal_number, synthesis_date_unixepoch,\
         initiator, monomer, terminator, initiator_long, monomer_long, terminator_long, long_name, product_A_smiles,\
          product_B_smiles, product_C_smiles, product_D_smiles, product_E_smiles, product_F_smiles, product_G_smiles,\
           product_H_smiles) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);',
        df[['vl_id_A', 'exp_nr', 'plate', 'well', 'lab_journal_nr', 'synthesis_date_unixepoch',
            'initiator', 'monomer', 'terminator', 'initiator_long', 'monomer_long', 'terminator_long', 'long_name',
            'SMILES_A',
            'SMILES_B', 'SMILES_C', 'SMILES_D', 'SMILES_E', 'SMILES_F', 'SMILES_G', 'SMILES_H']].values.tolist())
    con.commit()
    return


def main():
    con = sqlite3.connect(DB_PATH)

    # import the plate info
    plates_dict = get_plates_for_experiment(experiment_directory)

    # get the info about well <-> composition into a single df
    well_list = []
    for i, plate in plates_dict.items():
        for well, compounds, volume in plate.iterate_wells():
            well_list.append([i, well] + [c for c in compounds])
    well_df = pd.DataFrame(well_list, columns=['plate', 'well', 'initiator', 'monomer', 'terminator'])
    well_df = well_df.dropna()

    # get unique building blocks
    initiators = well_df['initiator'].unique().tolist()
    monomers = well_df['monomer'].unique().tolist()
    terminators = well_df['terminator'].unique().tolist()

    # get product infos from database
    print('Retrieving product data from virtuallibrary table...')
    products = bulk_data_retrieval_from_virtuallibrary(con, initiators, monomers, terminators)
    products = pd.DataFrame(products, columns=['vl_id',
                                               'initiator', 'monomer', 'terminator', 'initiator_long',
                                               'monomer_long', 'terminator_long', 'long_name', 'type', 'SMILES']) \
        .pivot(index=['initiator', 'monomer', 'terminator', 'initiator_long', 'monomer_long',
                      'terminator_long', 'long_name'], columns=['type', ],
               values=['vl_id', 'SMILES']) \
        .reset_index()  # pivoting is needed to go from a df where every type of product has its own row to a df where
    # one row holds all data for one well
    products.columns = ['_'.join(t).strip('_') for t in products.columns]  # pivoting leads to MultiIndex. Get rid of it

    # merge well and product info
    new_df = pd.merge(products, well_df, how='left', on=['initiator', 'monomer', 'terminator'])

    # add experiment info
    new_df['exp_nr'] = experiment_number
    new_df['lab_journal_nr'] = new_df['plate'].apply(lambda x: labj_nr_dict[x])
    new_df['synthesis_date_unixepoch'] = synthesis_date
    new_df = new_df.sort_values(['plate', 'well'])  # we sort so that DB table is a little more intuitive

    # write gathered info to experiments table
    print('Now writing to DB....')
    bulk_data_insertion_to_experiments(con, new_df)
    print('Finished writing to DB.')
    con.close()
    return


if __name__ == '__main__':
    main()
