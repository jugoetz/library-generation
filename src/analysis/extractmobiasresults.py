"""
Evaluate yields of a plate from MoBiAS output. Save the yields to the `lcms` table of the DB.

Requires the MoBiAS output file (`BMII00????_Skript-Results.csv`),
optionally output files from reprocessing (`BMII00????_Skript-Results_reprocess*.csv`) to overwrite selected records,
and the submission file (`BMIIyyyyyy-SampleTable_<lab journal number>.xls`) in the directory for the plate.

Give one or more lab journal numbers as an argument to the script (or in the config file). Only these will be processed.

This script will OVERWRITE existing records in the `lcms` table (with the same `synthesis_id`).
"""
import argparse
import re
import sqlite3

import pandas as pd
import numpy as np

from src.definitions import PLATES_DIR, DB_PATH
from src.util.utils import get_lcms_file_name, get_conf

# configuration
# edit config.yaml to change
conf = get_conf()


def import_lcms_results(path):
    """
    MoBiAS results come in a CSV file. Import the MoBiAS output into a dataframe, extract the sample ID,
    and ensure whitespaces in column names haven't been tampered with by accident.
    If we have additional result files from reprocessing (in a file named BMII00????_Skript-Results_reprocess*.csv),
    these records will overwrite the records from the original results file.
    """
    # read data
    df = pd.read_csv(
        path, header=3, encoding="latin-1", skip_blank_lines=False
    ).set_index(
        "Vial Pos"
    )  # read results.csv file from MoBiAS. Note that we use the Vial Pos index to match with reprocessed files if they exist
    # check if there are any reprocessed results files
    reprocessed_files = path.parent.glob(f"BMII00????_Skript-Results_reprocess*.csv")
    for reprocessed_file in reprocessed_files:
        df_reprocessed = pd.read_csv(
            reprocessed_file, header=3, encoding="latin-1", skip_blank_lines=False
        ).set_index("Vial Pos")
        # overwrite the original results with the reprocessed results
        df.update(df_reprocessed)
    # clean the sample ID. We only need everything after the last whitespace (the JG???-??? part)
    df["Sample ID"] = df["Sample ID"].str.split(" ").str[-1]
    # remove any accidental whitespaces in column names
    df.columns = df.columns.str.strip()
    return df.reset_index()  # put Vial Pos back as a column


def check_mobias_measurements_align_with_input(df, mobias_input, exp_nr):
    """
    Check that we received a measurement for all wells that were in the input and that the sum formulae are identical.
    (Basically a test against copy-paste or accidental editing errors)

    Print warnings for missing measurements and for sum formulae that differ between MoBiAS output and submission file.
    """

    # find relevant columns (SumFx) and discard everything else
    regex = re.compile("^SumF([0-9]+)$")

    area_columns = [column for column in df.columns if regex.match(column)]
    df.set_index("Sample ID", inplace=True)
    df = df[area_columns]

    # import the submission file and narrow it down to the same format as df
    df_input = pd.read_excel(mobias_input)
    df_input.set_index("Sample-Ident", inplace=True)
    for col in area_columns:
        if col not in df_input.columns:
            raise KeyError(f"Column {col} was not found in MoBiAS submission file")
    df_input = df_input[area_columns]
    df_input.drop("Blank", inplace=True)
    # df_input has empty values as float('nan'), df has them as '-'. Set all to '-'
    df_input.replace({float("nan"): "-"}, inplace=True)

    # check equality of df and df_input
    if np.all(df.eq(df_input)):
        # we are fine
        pass
    else:
        missing_indices = set([f"{exp_nr}-{i:03d}" for i in range(1, 321)]) - set(
            df.index.to_list()
        )
        if len(missing_indices) > 0:
            # some measurements are missing. We warn about this, but continue.
            print(
                f"Measurements missing from MoBiAS output. These measurements are missing in the "
                f"output: {missing_indices}"
            )

        else:
            # no measurements missing, but some sum formulae are different. (This may happen due to reprocessing)
            # check that the sum formulae are the same, and warn if not
            for index in df.index:
                if not np.all(df.loc[index].eq(df_input.loc[index])):
                    print(
                        f"Warning: Sum formulae for {index} differ between MoBiAS output and submission file. "
                        f"MoBiAS: {df.loc[index].to_list()}, Submission: {df_input.loc[index].to_list()}"
                    )

    return


def clean_result_df(df):
    """Discard the unneeded columns, adjust dtypes and extract row and column of the analyzed wells"""

    # Remove unneeded columns. We want to keep Vial Pos, File (the raw data name) and all MS Areas
    columns = ["Vial Pos", "File"]
    columns += [
        s for s in df.columns if "Area" in s and "UV" not in s and df[s].count() > 0
    ]
    df = df.loc[:, columns]

    # For fields where we did not search for a compound, '-' appears. Replace with np.nan.
    df.replace(to_replace="-", value=np.nan, inplace=True)
    # Since '-' would have induced object dtype, we change to float64. We do this for all Area columns, to be safe.
    area_columns = [c for c in df.columns if c.endswith("Area")]
    df[area_columns] = df[area_columns].astype("float64")
    # split "Vial Pos" into separate columns for plate, row, column. Discard "Vial Pos"
    df["plate"] = df["Vial Pos"].str.split("-").str[0]
    df["row"] = df["Vial Pos"].str.split("-").str[1]
    df["column"] = df["Vial Pos"].str.split("-").str[2]
    df.drop(columns=["Vial Pos"], inplace=True)

    return df


def save_mobias_data_to_db(df, db_path, exp_nr):
    """
    Create a table 'lcms' to hold raw lcms results (if it does not exist). Save the results from dataframe to DB.
    The data is split into one column for identities (SumF1 etc.) and one for areas (e.g. 36574.000).
    Both columns contain str-representations of lists. compounds[0] corresponds to area[0]
    """
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA foreign_keys = 1")
    cur = con.cursor()
    # if the table for LCMS data does not exist, create it
    cur.execute(
        "CREATE TABLE IF NOT EXISTS lcms ("
        "id INTEGER PRIMARY KEY ,"
        "synthesis_id INTEGER,"
        "lcms_compounds TEXT,"
        "lcms_areas TEXT,"
        "FOREIGN KEY(synthesis_id) REFERENCES experiments(id)"
        ");"
    )
    con.commit()
    for i, row in df.iterrows():
        well = f'{row["row"]}{row["column"]}'
        try:
            synthesis_id = cur.execute(
                "SELECT id FROM experiments WHERE lab_journal_number = ? AND well = ?;",
                (exp_nr, well),
            ).fetchone()[0]
        except TypeError:
            raise LookupError(
                f'No entry in database table "experiments" for lab journal number {exp_nr} and well {well}'
            )
        compounds, areas = [], []
        for index, value in row.items():
            if index.endswith("Area"):
                compounds.append(index)
                areas.append(value)
        cur.execute(
            "INSERT INTO lcms (synthesis_id, lcms_compounds, lcms_areas) VALUES (?, ?, ?) ON CONFLICT(synthesis_id) DO UPDATE SET lcms_compounds = ?, lcms_areas = ? WHERE synthesis_id = ?;",
            (
                synthesis_id,
                repr(compounds),
                repr(areas),
                repr(compounds),
                repr(areas),
                synthesis_id,
            ),
        )
    con.commit()
    return


def extract_mobias_results(path, db_path, exp_nr, mobias_input):
    """
    Main function. Extract the csv supplied by MoBiAS, verify against my submission file, clean up the data,
    and save the extracted data to a database.
    """
    results_df = import_lcms_results(path)
    check_mobias_measurements_align_with_input(results_df, mobias_input, exp_nr)
    results_df = clean_result_df(results_df)
    save_mobias_data_to_db(results_df, db_path, exp_nr)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "lab_journal_numbers",
        type=str,
        default=conf["lab_journal_numbers"],
        nargs="*",
        help="Lab journal numbers to process",
    )
    args = parser.parse_args()

    for lcms_nr in args.lab_journal_numbers:
        exp_dir = PLATES_DIR / lcms_nr
        results_file_name = get_lcms_file_name(lcms_nr)

        print(f"Now extracting data for {lcms_nr} from {results_file_name}...")

        results_file_path = exp_dir / results_file_name
        submission_file = exp_dir / f"BMIIyyyyyy-SampleTable_{lcms_nr}.xls"
        extract_mobias_results(results_file_path, DB_PATH, lcms_nr, submission_file)
    print("Finished data extraction!")
