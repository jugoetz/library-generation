"""
Evaluate yields of a plate from Mobias output. Save the yields to the DB.
"""

import pandas as pd
from config import *
import sqlite3
import numpy as np
import re

"""GLOBALS"""
# TODO having these is not ideal. Better to read them from DB or infer them from os.walk
results_file_path = EXP_DIR / 'BMII001987_Skript-Results_protecting-groups.csv'
submission_file = EXP_DIR / 'BMIIyyyyyy-SampleTable_JG216-big.xls'
with open(EXP_DIR / 'notebook_nr.txt') as file:
    LCMS_number = file.read().strip('\n').strip()


def import_lcms_results(path):
    """
    Mobias results come in a CSV file. Import the Mobias output into a dataframe, extract the sample ID,
    and ensure Mobias boys haven't been messing with us again by adding whitespaces to the column names.
    """
    # read data
    df = pd.read_csv(path, header=3, encoding='latin-1')  # read results.csv file from Mobias
    # extract the sample ID (the JG2xx-001 part). It is part after the last whitespace in that field
    df["Sample ID"] = df["Sample ID"].str.split(" ").str[-1]
    # remove any accidental whitespaces in column names
    df.columns = df.columns.str.strip()
    return df


def check_mobias_input_output_equivalent(df, mobias_input):
    """Check if the sum formulae in the Mobias output are the same that I had entered in my submission."""

    # find relevant columns (SumFx) and discard everything else
    regex = re.compile('^SumF([0-9]+)$')

    area_columns = [column for column in df.columns if regex.match(column)]
    df.set_index('Sample ID', inplace=True)
    df = df[area_columns]

    # import the submission file and narrow it down to the same format as df
    df_input = pd.read_excel(mobias_input)
    df_input.set_index('Sample-Ident', inplace=True)
    for col in area_columns:
        if col not in df_input.columns:
            raise KeyError(f'Column {col} was not found in mobias submission file')
    df_input = df_input[area_columns]
    df_input.drop('Blank', inplace=True)
    # df_input has empty values as float('nan'), df has them as '-'. Set all to '-'
    df_input.replace({float('nan'): '-'}, inplace=True)

    # check equality of df and df_input
    if np.alltrue(df.eq(df_input)):
        # we are fine
        pass
    else:
        raise ValueError('Values for Mobias input and output do not align')

    return


def clean_result_df(df):
    """Discard the unneeded columns, adjust dtypes and extract row and column of the analyzed wells"""

    # Remove unneeded columns. We want to keep Vial Pos, File (the raw data name) and all MS Areas
    columns = ['Vial Pos', 'File']
    columns += [s for s in df.columns if "Area" in s and "UV" not in s and df[s].count() > 0]
    df = df[columns]

    # For fields where we did not search for a compound, '-' appears. Replace with np.nan.
    df.replace(to_replace='-', value=np.nan, inplace=True)
    # Since '-' would have induced object dtype, we change to float64. We do this for all Area columns, to be save.
    for c in df.columns:
        if c.endswith('Area'):
            df[c] = df[c].astype('float64')

    # split "Vial Pos" into separate columns for plate, row, column. Discard "Vial Pos"
    df['plate'] = df['Vial Pos'].str.split('-').str[0]
    df['row'] = df['Vial Pos'].str.split('-').str[1]
    df['column'] = df['Vial Pos'].str.split('-').str[2]
    df.drop(columns=['Vial Pos'], inplace=True)

    return df


def save_mobias_data_to_db(df, db_path, exp_nr):
    """
    Create a table 'lcms' to hold raw lcms results (if it does not exist). Save the results from dataframe to DB.
    The data is split into one column for identities (SumF1 etc.) and one for areas (e.g. 36574.000).
    Both columns contain str-representations of lists. compounds[0] corresponds to area[0]
    """
    con = sqlite3.connect(db_path)
    con.execute('PRAGMA foreign_keys = 1')
    cur = con.cursor()
    # if the table for LCMS data does not exist, create it
    cur.execute('CREATE TABLE IF NOT EXISTS lcms ('
                'id INTEGER PRIMARY KEY ,'
                'synthesis_id INTEGER,'
                'lcms_compounds TEXT,'
                'lcms_areas TEXT,'
                'FOREIGN KEY(synthesis_id) REFERENCES experiments(id)'
                ');')
    con.commit()
    for i, row in df.iterrows():
        well = f'{row["row"]}{row["column"]}'
        synthesis_id = \
            cur.execute('SELECT id FROM experiments WHERE lab_journal_number = ? AND well = ?;',
                        (exp_nr, well)).fetchone()[
                0]
        compounds, areas = [], []
        for index, value in row.iteritems():
            if index.endswith('Area'):
                compounds.append(index)
                areas.append(value)
        cur.execute('INSERT INTO lcms (synthesis_id, lcms_compounds, lcms_areas) VALUES (?, ?, ?);',
                    (synthesis_id, repr(compounds), repr(areas)))
    con.commit()
    return


def extract_mobias_results(path, db_path, exp_nr, mobias_input):
    """
    Main function. Extract the csv supplied by mobias, verify against my submission file, clean up the data,
    and save the extracted data to a database.
    """
    results_df = import_lcms_results(path)
    check_mobias_input_output_equivalent(results_df, mobias_input)
    results_df = clean_result_df(results_df)
    save_mobias_data_to_db(results_df, db_path, exp_nr)
    return


if __name__ == '__main__':
    extract_mobias_results(results_file_path, DB_PATH, EXP_NR, submission_file)