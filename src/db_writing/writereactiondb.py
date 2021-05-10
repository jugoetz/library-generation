import sqlite3
from config import *
from datetime import datetime
from db_retrieval.generatelcmssubmission import import_pl
import os

con = sqlite3.connect(DB_PATH)
cur = con.cursor()

"""Some functions for simple queries to the db"""


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
    # control variables
    experiment_number = 2
    synthesis_date = datetime(2021, 5, 11).timestamp()

    plate_nr_to_labj_nr = {
        '1': 'JG222',
        '2': 'JG223',
        '3': 'JG224',
        '4': 'JG225',
        '5': 'JG226',
        '6': 'JG227',
    }
    plates_dict = {}
    for path, _, files in os.walk(PLATES_DIR / f'exp{experiment_number}'):
        for f in files:
            m = PLATE_REGEX.match(f)
            if m:
                plate_dict = import_pl(Path(path, f))
                plates_dict[m.group(1)] = plate_dict

    for plate_number, plate_content in plates_dict.items():
        labj_nr = plate_nr_to_labj_nr[plate_number]
        for well, content in plate_content.items():
            if type(content) is list and len(content) == 3:
                copy_data(content[0], content[1], content[2], experiment_number, plate_number, well, labj_nr,
                          synthesis_date)
    con.close()
