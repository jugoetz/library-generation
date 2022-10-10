"""
Populate the 'experiments' table of the database.
Location data is taken from plate_layout files in either experiment directories (e.g. exp2, standard usage) or
in individual directories with lab journal numbers (e.g. JG228, non-standard - add things not within 50k project scope).
Compound data is retrieved from 'virtuallibrary' table.

Note that this script appends records to the 'experiments' table.
Any previous, outdated records will remain in the table.
"""

import re
import os
import csv
import sqlite3
import shutil

import pandas as pd

from src.definitions import PLATES_DIR, DB_PATH, PLATE_LIST_PATH
from src.experiment_planning.generatelcmssubmission import import_pl
from src.util.utils import get_conf

# configuration
# edit config.yaml to change
conf = get_conf()
if len(conf["db"]["plate_nrs"]) == len(conf["lab_journal_numbers"]):
    conf["db"]["lab_journal_nr_dict"] = dict(
        zip(conf["db"]["plate_nrs"], conf["lab_journal_numbers"])
    )
else:
    raise ValueError(
        "db.plate_nrs and lab_journal_numbers do not align (different length)"
    )


def make_new_directories(dir_list):
    """Make a list of directories if they don't exist yet and copy the respective plate layout to them"""
    print(f"Making new directories (if they do not exist): {repr(dir_list)}")
    for i in dir_list:
        (PLATES_DIR / i).mkdir(exist_ok=True)
    return


def copy_plate_layout_files(exp_dir, mapping_dict):
    """Copy the plate layout files from the experiment directory (e.g.exp1) to the plate directories (e.g JG255)"""
    for src_nr, target_dir in mapping_dict.items():
        src_main = PLATES_DIR / exp_dir / f"plate_layout_plate{src_nr}.csv"
        src_vol = PLATES_DIR / exp_dir / f"plate_layout_plate{src_nr}_volumes.csv"
        target = PLATES_DIR / target_dir
        print(f"Copying \t{src_main}\t to \t{target}")
        print(f"Copying \t{src_vol}\t to \t{target}")
        shutil.copy2(src_main, target)
        shutil.copy2(src_vol, target)
    return


def append_plate_log(exp_dir, exp_nr, mapping_dict):
    """Append the log in plate_list.csv with mapping information for the current plates"""
    with open(PLATE_LIST_PATH, "a") as file:
        writer = csv.writer(file)
        writer.writerows(
            [
                [exp_dir, exp_nr, plate_nr, labj_nr, ""]
                for plate_nr, labj_nr in mapping_dict.items()
            ]
        )
    return


def get_plates_for_experiment(exp_dir):
    """
    Import the plate layouts for an experiment into a dictionary of plates.
    """
    plates_dict = {}
    with os.scandir(PLATES_DIR / f"{exp_dir}") as it:
        for entry in it:
            m = re.compile(conf["plate_regex"]).match(entry.name)
            if m:
                print(f"Retrieving plate layout: {entry.name}...")
                plate = import_pl(entry.path, return_type="plate")
                plates_dict[m.group(1)] = plate
    print(f"Retrieved data for {len(plates_dict)} plate(s).")
    return plates_dict


def bulk_data_retrieval_from_virtuallibrary(
    con, initiators: list, monomers: list, terminators: list
):
    """
    From the virtuallibrary table, retrieve all combinations of the all initiators,
    monomers and terminators given to the function.
    """
    cur = con.cursor()
    print("Retrieving product data from virtuallibrary table...")
    results = cur.execute(
        f"""SELECT id, initiator, monomer, terminator, initiator_long, monomer_long, terminator_long,
     long_name, type, SMILES FROM main.virtuallibrary
     WHERE initiator IN ({(", ".join("?" for _ in initiators))})
     AND monomer IN ({(", ".join("?" for _ in monomers))})
     AND terminator IN ({(", ".join("?" for _ in terminators))});""",
        initiators + monomers + terminators,
    ).fetchall()
    print(f"Retrieved product data for {len(results)} products.")
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
    print("Now writing reactions to DB....")
    cur.executemany(
        "INSERT INTO main.experiments (vl_id, exp_nr, plate_nr, well, lab_journal_number, synthesis_date_unixepoch,\
         initiator, monomer, terminator, initiator_long, monomer_long, terminator_long, long_name, product_A_smiles,\
          product_B_smiles, product_C_smiles, product_D_smiles, product_E_smiles, product_F_smiles, product_G_smiles,\
           product_H_smiles) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);",
        df[
            [
                "vl_id_A",
                "exp_nr",
                "plate",
                "well",
                "lab_journal_nr",
                "synthesis_date_unixepoch",
                "initiator",
                "monomer",
                "terminator",
                "initiator_long",
                "monomer_long",
                "terminator_long",
                "long_name",
                "SMILES_A",
                "SMILES_B",
                "SMILES_C",
                "SMILES_D",
                "SMILES_E",
                "SMILES_F",
                "SMILES_G",
                "SMILES_H",
            ]
        ].values.tolist(),
    )
    con.commit()
    print(f"Wrote {len(df)} reactions to DB.")
    return


def main():
    con = sqlite3.connect(DB_PATH)

    # make new experiment directories
    make_new_directories(list(conf["db"]["lab_journal_nr_dict"].values()))

    # copy plate files
    copy_plate_layout_files(conf["exp_dir"], conf["db"]["lab_journal_nr_dict"])

    # add the plate mapping information to plate_list.csv
    append_plate_log(conf["exp_dir"], conf["exp_nr"], conf["db"]["lab_journal_nr_dict"])

    # import the plate info
    plates_dict = get_plates_for_experiment(conf["exp_dir"])

    # we only consider plates listed in the config file
    plates_dict = {
        k: v for k, v in plates_dict.items() if int(k) in conf["db"]["plate_nrs"]
    }
    print(
        f"Based on the configuration for plate_nrs, {len(plates_dict)} plate(s) are considered."
    )

    # get the info about well <-> composition into a single df
    well_list = []
    for i, plate in plates_dict.items():
        for well, compounds, volume in plate.iterate_wells():
            well_list.append([i, well] + [c for c in compounds])
    well_df = pd.DataFrame(
        well_list, columns=["plate", "well", "initiator", "monomer", "terminator"]
    )
    well_df = well_df.dropna()

    # get unique building blocks
    initiators = well_df["initiator"].unique().tolist()
    monomers = well_df["monomer"].unique().tolist()
    terminators = well_df["terminator"].unique().tolist()

    # get product infos from database
    products = bulk_data_retrieval_from_virtuallibrary(
        con, initiators, monomers, terminators
    )
    products = (
        pd.DataFrame(
            products,
            columns=[
                "vl_id",
                "initiator",
                "monomer",
                "terminator",
                "initiator_long",
                "monomer_long",
                "terminator_long",
                "long_name",
                "type",
                "SMILES",
            ],
        )
        .pivot(
            index=[
                "initiator",
                "monomer",
                "terminator",
                "initiator_long",
                "monomer_long",
                "terminator_long",
                "long_name",
            ],
            columns=[
                "type",
            ],
            values=["vl_id", "SMILES"],
        )
        .reset_index()
    )  # pivoting is needed to go from a df where every type of product has its own row to a df where
    # one row holds all data for one well
    products.columns = [
        "_".join(t).strip("_") for t in products.columns
    ]  # pivoting leads to MultiIndex. Get rid of it

    print(
        f"After processing, products for {len(products)} unique reagent combinations have been obtained."
    )
    # merge well and product info
    # Note that this forms the intersection, so if entries are missing in one of the dfs, this might fail silently.
    # This behavior is desirable because sometimes, we do not exhaustively synthesise all combinations in 'products'.
    # To alert the user to silent failing, we print the number of reactions below.
    new_df = pd.merge(
        products, well_df, how="inner", on=["initiator", "monomer", "terminator"]
    )

    # add experiment info
    new_df["exp_nr"] = conf["exp_nr"]
    new_df["lab_journal_nr"] = new_df["plate"].apply(
        lambda x: conf["db"]["lab_journal_nr_dict"][int(x)]
    )
    new_df["synthesis_date_unixepoch"] = conf["db"]["synthesis_date"]
    new_df = new_df.sort_values(
        ["plate", "well"]
    )  # we sort so that DB table is a little more intuitive

    # write gathered info to experiments table
    bulk_data_insertion_to_experiments(con, new_df)
    con.close()
    return


if __name__ == "__main__":
    main()
