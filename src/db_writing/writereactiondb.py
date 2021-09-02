"""
Populate the 'experiments' table of the database.
Location data is taken from plate_layout files in either experiment directories (e.g. exp2, standard usage) or
in individual directories with lab journal numbers (e.g. JG228, non-standard - add things not within 50k project scope).
Compound data is retrieved from 'virtuallibrary' table.
"""

import re
import os
import sqlite3
from datetime import datetime
from pathlib import Path

from definitions import PLATES_DIR, DB_PATH
from db_retrieval.generatelcmssubmission import import_pl
from utils import get_conf

conf = get_conf()
con = sqlite3.connect(DB_PATH)
cur = con.cursor()

# control variables
experiment_directory = 'exp7'
experiment_number = 7  # can be int or None if not a canonical (50 k) plate
synthesis_date = datetime(2021, 9, 14).timestamp()
plate_nr_to_labj_nr = {
    '1': 'JG264',
    '2': 'JG265',
    '3': 'JG266',
    '4': 'JG267',
    '5': 'JG268',
    '6': 'JG269',
}


def copy_data(i, m, t, exp_nr, plate_nr, well, lab_journal_number, synthesis_date_unixepoch):
    # gather information from virtual_library database
    result = cur.execute(
        'SELECT id, initiator, monomer, terminator, initiator_long, monomer_long, terminator_long, long_name, type, SMILES FROM main.virtuallibrary WHERE initiator = ? AND monomer = ? AND terminator = ?;',
        (i, m, t)).fetchall()
    vl_id, initiator, monomer, terminator, initiator_long, monomer_long, terminator_long, long_name, _, _ = result[0]
    # rearrange smiles of the different query hits
    smiles = {i[-2]: i[-1] for i in result}  # gives dict of <product type>: <smiles>
    # if any product is unassigned, assign a string NA
    for t in 'ABCDEFGH':
        if t not in smiles.keys():
            smiles[t] = 'NA'
    # write gathered information to experiments database
    cur.execute(
        'INSERT INTO main.experiments (vl_id, exp_nr, plate_nr, well, lab_journal_number, synthesis_date_unixepoch, initiator, monomer, terminator, initiator_long, monomer_long, terminator_long, long_name, product_A_smiles, product_B_smiles, product_C_smiles, product_D_smiles, product_E_smiles, product_F_smiles, product_G_smiles, product_H_smiles) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);',
        (vl_id, exp_nr, plate_nr, well, lab_journal_number, synthesis_date_unixepoch, initiator, monomer, terminator,
         initiator_long, monomer_long, terminator_long, long_name, smiles['A'], smiles['B'], smiles['C'], smiles['D'],
         smiles['E'], smiles['F'], smiles['G'], smiles['H']))
    con.commit()
    return


if __name__ == '__main__':
    # import the plate layouts, i.e. retrieve well location data
    plates_dict = {}
    for path, _, files in os.walk(PLATES_DIR / f'{experiment_directory}'):
        for f in files:
            m = re.compile(conf['plate_regex']).match(f)
            if m:
                print(f'Retrieving plate layout {f}...')
                plate_dict = import_pl(Path(path, f))
                plates_dict[m.group(1)] = plate_dict

    # retrieve compound data from 'virtuallibrary' table and write to 'experiments' table
    # TODO these two steps should be separated(?) and the call to DB probably could use optimization
    for plate_number, plate_content in plates_dict.items():
        print(f'Retrieving reaction data for plate {plate_number} and saving to DB...')
        labj_nr = plate_nr_to_labj_nr[plate_number]
        for well, content in plate_content.items():
            if type(content) is list and len(content) == 3:
                copy_data(content[0], content[1], content[2], experiment_number, plate_number, well, labj_nr,
                          synthesis_date)
    con.close()
