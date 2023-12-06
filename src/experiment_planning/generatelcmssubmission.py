"""
Generate the submission file for high-throughput MS of library compounds.

Inputs:
  - Plate layout (csv-files, usually 'plate_layout_plateX.csv'):
    One plate per file. Handles arbitrary number of plates/files. X may be any number.
  - EITHER
     the mol_prop_dict prepared from sdf_to_properties (legacy)
    OR
     a SQLite DB table of the virtual library (standard now)

  If shorthand names are used in the plate layout, the long names are read from the DB.

Output:
  - Submission file (mobias_submission.csv):
    One well per row.
    Columns are:
      1 well identifier (e.g. Py-A-3)
    + 8 molecular formulae for SynFerm products A-H
    + a variable number of deprotection product columns (e.g. from Boc removal)
    + 1 column for internal standard
"""
import argparse
import json
import re
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd

from src.definitions import PLATES_DIR
from labware.plates import Plate384, Plate96
from src.util.utils import get_conf
from src.util.db_utils import SynFermDatabaseConnection

# configuration
# edit config.yaml to change
conf = get_conf()
debug = False
con = SynFermDatabaseConnection()


def import_pl(file, return_type="dict"):
    """
    Import the plate layout from x files 'plate_layout_plateX.csv' where X is the plate number.
    Return dict of the form {'A1': ['I1', 'M1', 'T1']} that maps wells to shorthand names if return_type == 'dict'
    OR
    Return labware.plates.Plate if return_type == 'plate'
    :return: dict or labware.plates.Plate
    """
    if conf["lcms"]["well_plate_size"] == 96:
        p = Plate96(150000, 5000)
    elif conf["lcms"]["well_plate_size"] == 384:
        p = Plate384(12000, 2500)
    else:
        raise ValueError(f'Invalid plate size: {conf["lcms"]["well_plate_size"]}')
    p.from_csv(file)
    if return_type == "dict":
        return p.to_dict()
    elif return_type == "plate":
        return p
    else:
        raise ValueError(
            f'Unknown return_type {return_type}. return_type can be "plate" or "dict".'
        )


def get_long_name(row, exp_nr):
    """
    For use with pandas.DataFrame.apply()
    Convert shorthand names located in different df columns (I, M, T) to one longhand name.
    """
    long = []
    for col in ["I", "M", "T"]:
        if (
            row[col] is None
            or row[col] == "None"
            or row[col] == ""
            or row[col] == "n/a"
        ):
            pass  # catch empty fields (this will usually happen)
        else:
            if row[col][1:].isdigit():  # if we have short names, get the long name
                # if the column contains a long name already, just append it
                con.get_long_name(short=row[col], exp_nr=exp_nr)
            else:  # if we have long names, we use that unchanged
                long.append(row[col])
    return " + ".join(long)


def get_prop_from_db():
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
    :return: dict
    """
    cur = con.con.cursor()
    prop_dict = {}
    for row in cur.execute(
        "SELECT id, long_name, type, molecular_formula_1, molecular_formula_alt, lcms_mass_1, lcms_mass_alt FROM virtuallibrary"
    ).fetchall():
        if (
            not row[1] in prop_dict
        ):  # if we have not encountered this long name before, make a dictionary as the entry
            prop_dict[row[1]] = {}
        prop_dict[row[1]][row[2]] = (
            row[3],
            row[5],
        )  # these are the "standard" formula and mass
        if row[4] is not None:  # if alternate masses are present
            # handle different input formats
            formulas = [
                s.strip("[' ]") for s in row[4].split(",")
            ]  # works for CH4 and ['CH4'] and ['CH4', 'C2H6']
            if isinstance(row[6], str):
                masses = [
                    float(s.strip().strip("[' ]")) for s in row[6].split(",")
                ]  # works for "150.0, 10.0" and ['150.0', '10.0']
            elif isinstance(row[6], float):
                masses = [row[6]]  # works for 150.0
            else:
                raise RuntimeError(f"Received unexpected type {type(row[4])}")
            for i, (f, m) in enumerate(zip(formulas, masses)):
                prop_dict[row[1]][f"{row[2]}_{i + 2}"] = (f, m)
    return prop_dict


def get_prop(long_name, mol_props, prop):
    """
    Get mass and formula of a molecule identified by long_name from the mol_props dictionary.
    :param long_name: str
    :param mol_props: dict
    :param prop: str, either of "mass" or "formula"
    """
    if prop == "mass":
        i = 0
    elif prop == "formula":
        i = 1
    else:
        raise ValueError("")
    if long_name == "":
        return None
    else:
        try:
            mass = mol_props[long_name][i]
        except KeyError:
            mass = "n/a"
    return mass


def add_is(df):
    df["IS_mass"] = np.nan
    df["IS_formula"] = ""
    df.loc[df["long"] != "", "IS_mass"] = conf["lcms"]["is_mass"]
    df.loc[df["long"] != "", "IS_formula"] = conf["lcms"]["is_formula"]
    return df


def write_csv(df, file, exp_dir):
    """
    Generate formatted csv output for MoBiAS
    """

    def splitwell(df):
        well = df["well"]
        new = f"Py-{str(well)[0]}-{str(well)[1:]}"
        return new

    if conf["lcms"]["mass_or_formula"] == "mass":
        suffix = "_mass"
    elif conf["lcms"]["mass_or_formula"] == "formula":
        suffix = "_formula"
    else:
        raise ValueError(f'Invalid option {conf["lcms"]["mass_or_formula"]}')
    df["Vial"] = df.loc[:, ["plate", "well"]].apply(splitwell, axis=1)
    # drop any row were all masses are np.nan
    df.dropna(
        axis=0,
        how="all",
        subset=[col for col in df.columns if col.endswith("mass") and col != "IS_mass"],
        inplace=True,
    )
    subset = [
        "Vial",
    ]

    def sort_key(s):
        if s == "Vial":
            primary = 0
        elif s.startswith("IS"):
            primary = 3
        elif any(c.isdigit() for c in s):
            primary = 2
        else:
            primary = 1
        return (primary, s)

    for s in df.columns:
        if s.endswith(suffix):
            subset.append(s)
    # we order the subset list according to: Vial, A-H, A2-Ax...H2-Hx, IS
    subset = sorted(subset, key=sort_key)
    subset_df = df.loc[
        df["long"] != "", subset
    ]  # remove columns with no long name set (those are emtpy wells e.g. A1)
    # define rename dict for column renaming.
    # should always map product A_x -> SumF1 ... H_x -> SumF8, IS -> SumF9
    rename_dict = {}
    counter = 1
    for column in subset_df.columns:
        if column.endswith(suffix):
            rename_dict[column] = f"SumF{counter}"
            counter += 1

    subset_df.rename(columns=rename_dict, inplace=True)
    with open(exp_dir / "compound_alternative_mass_dict.json", "w") as jsonfile:
        json.dump(rename_dict, jsonfile)  # save this

    subset_df.replace(
        "n/a", "", inplace=True
    )  # we have n/a strings where enumeration has not given a product. Replace them with empty string.
    # finally, resort the dataframe columns so that vial is in front and SumFx ascending
    cols = subset_df.columns.tolist()
    cols.remove("Vial")
    reordered_columns = ["Vial"] + sorted(cols, key=lambda x: int(x.strip("SumF")))
    subset_df = subset_df[reordered_columns]
    subset_df.to_csv(file, index=False)
    return


def main(exp_dir):
    # Import Mass and Formula from DB
    mol_prop_dict = get_prop_from_db()

    # Import plates
    plates_dict = {}
    for path, _, files in os.walk(exp_dir):
        for f in files:
            m = re.compile(conf["plate_regex"]).match(f)
            if m:
                plate_dict = import_pl(Path(path, f))
                plates_dict[m.group(1)] = plate_dict

    # Put everything in df for ease of use
    dfs = []
    for plate, plate_content in plates_dict.items():
        df = pd.DataFrame.from_dict(
            plate_content, orient="index", columns=["I", "M", "T"]
        )
        df.reset_index(drop=False, inplace=True)
        df.rename({"index": "well"}, axis=1, inplace=True)
        df.insert(loc=0, column="plate", value=int(plate))
        dfs.append(df)
    df = pd.concat(dfs)
    df.sort_values(by=["plate"], inplace=True, kind="mergesort")  # mergesort
    df.reset_index(drop=True, inplace=True)

    # Translate product names from shorthand to long names (e.g. 'Al002 + Mon001 + TerTH010')
    df["long"] = df.loc[:, ["I", "M", "T"]].apply(
        get_long_name, axis=1, exp_nr=con.get_exp_nr(exp_dir.name)
    )

    if debug:
        # print the input values for double-checking
        print("########## INPUT VALUES ###########\n")
        for k, v in mol_prop_dict.items():
            print(f"Products {k}:\n{v}\n")
        print(
            f"Used starting materials: \n{np.unique(df.loc[:, ['I', 'M', 'T']].to_numpy().flatten()).tolist()}\n\n"
        )
        for k, v in plates_dict.items():
            print(f"Plate layout {k}:\n{v}\n")

        # start printing terminal output
        print("########## OUTPUT ###########\n")
        for i, data in df.iterrows():
            print(
                f'Plate {data["plate"]}, Well {data["well"]}, Product: {data["long"]}'
            )

    # Add molecular formulae and exact masses to dataframe:
    # select the subportion of the prop-dict that we need (for this specific plate)
    needed_dict = {k: v for k, v in mol_prop_dict.items() if k in df["long"].values}
    # find how many different submission forms we will need
    # Here's how this works: We only need look at products A and E since they are candidates for the
    # maximum number of combinations. Everything else will have equal or less.
    # We go through the keys in the inner dict-level, only take unique ones starting with A, resp. D. The length
    # of this set corresponds to the maximum number of combinations and thus submission forms

    combinations_a = len(
        set(
            [
                key
                for v in needed_dict.values()
                for key in v.keys()
                if key.startswith("A")
            ]
        )
    )
    combinations_d = len(
        set(
            [
                key
                for v in needed_dict.values()
                for key in v.keys()
                if key.startswith("D")
            ]
        )
    )
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
                if index == 1 and letter in ["A", "B", "C", "D", "E", "F", "G", "H"]:
                    # if the columns are not yet present, fill them with placeholder values
                    if f"{letter}_mass" not in df_this.columns:
                        df_this[f"{letter}_mass"] = np.nan
                    if f"{letter}_formula" not in df_this.columns:
                        df_this[f"{letter}_formula"] = ""
                    df_this.loc[
                        df_this["long"] == long_name, f"{letter}_mass"
                    ] = mol_props[1]
                    df_this.loc[
                        df_this["long"] == long_name, f"{letter}_formula"
                    ] = mol_props[0]
                if letter.endswith(str(index)):
                    # if the columns are not yet present, fill them with placeholder values
                    if f"{letter}_mass" not in df_this.columns:
                        df_this[f"{letter}_mass"] = np.nan
                    if f"{letter}_formula" not in df_this.columns:
                        df_this[f"{letter}_formula"] = ""
                    df_this.loc[
                        df_this["long"] == long_name, f"{letter}_mass"
                    ] = mol_props[1]
                    df_this.loc[
                        df_this["long"] == long_name, f"{letter}_formula"
                    ] = mol_props[0]

    # from the individual dataframes for different -PG product sets, we form one big dataframe and drop the na columns
    df = dfs[0]
    for i in dfs[1:]:
        df = pd.merge(df, i, how="left")

    # add internal standard if user wishes
    if conf["lcms"]["add_is"] is True:
        add_is(df)
    else:
        print("You chose not to add internal standard.\n")

    # save to file
    output_file = exp_dir / "mobias_submission.csv"
    write_csv(df, output_file, exp_dir)
    print(f'Data was written to "{output_file}".')

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

    for lab_journal_nr in args.lab_journal_numbers:
        print(f"Generating submission file for {lab_journal_nr}...")
        main(PLATES_DIR / lab_journal_nr)
    print("End of script. Exiting...")
