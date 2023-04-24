"""
This script is supposed to run after all experiment data has been collected.
It should give a list of plates/wells that had errors and save it to the experiment directory.
Further it should save this information to the experiments DB table.
A column "valid" is supposed to hold the information. Its value will be NULL for valid experiments or a verbose
description of the error for potentially invalid experiments.
Verbose descriptions contain at least one of the keywords WARNING or ERROR depending on the severity of the issue.

The following are considered ERRORS:
- A manually curated list of errors, e.g. to record pipetting errors / mix-ups / precipitation etc.
- Transfer errors at NEXUS extracted from the provided transfer files
- Wells with low volume (< 2.2 uL) in last NEXUS survey after overnight incubation
- More than 4 peaks for the main product in MoBiAS analysis
- More than 1 peak for the internal standard in MoBiAS analysis

The following are considered WARNINGS::
- 2-4 peaks for the main product in MoBiAS analysis
- IS peak area <50% of the median for the plate
- IS peak area >200% of the median for the plate

Note that the script will overwrite the contents of the "valid" column in the DB table.
"""
import argparse
import sqlite3
import warnings
from pathlib import Path
from typing import List

import pandas as pd

from src.definitions import PLATES_DIR, PLATE_LIST_PATH, DB_PATH
from src.util.utils import get_internal_standard_number, get_conf

# configuration
# edit config.yaml to change
conf = get_conf()


def get_manual_error_records(path: Path) -> List[str]:
    """
    Import a CSV file containing manually curated error records.
    :param path: path to manual error file (csv)
    :return: List of manual errors
    """
    try:
        errors = pd.read_csv(path, dtype=object)
    except FileNotFoundError:
        warnings.warn(
            f"No manual error list found in {path.parent}. Manual errors will not be included."
        )
        return []
    return errors.to_numpy(dtype="str").tolist()


def read_nexus_transfer_errors(path: Path, error_string: str) -> List[str]:
    """
    Import a CSV transfer file and add everything under exceptions to error records.
    :param path: path to transfer file (csv)
    :param error_string: verbose description of the error to add to the error list
    :return: List of transfer errors
    """
    transfers = pd.read_csv(
        path,
        header=None,
        names=[
            "Source Plate Name",
            "Source Plate Barcode",
            "Source Plate Type",
            "Source Well",
            "Source Concentration",
            "Source Concentration Units",
            "Destination Plate Name",
            "Destination Plate Barcode",
            "Destination Well",
            "Destination Concentration",
            "Destination Concentration Units",
            "Compound Name",
            "Transfer Volume",
            "Actual Volume",
            "Transfer Status",
            "Current Fluid Height",
            "Current Fluid Volume",
            "% DMSO",
        ],
    )
    if "[EXCEPTIONS]" in transfers.values:
        exceptions_idx = transfers.loc[
            transfers["Source Plate Name"] == "[EXCEPTIONS]"
        ].index
        details_idx = transfers.loc[transfers["Source Plate Name"] == "[DETAILS]"].index
        transfers["plate"] = (
            transfers["Destination Plate Barcode"]
            .str.strip("Synthesis")
            .str.strip("Analysis")
        )
        # ^ we can do this since the analysis plate is a 1:1 copy
        transfers["row"] = transfers["Destination Well"].str[0]
        transfers["column"] = transfers["Destination Well"].str[1:]
        transfers["error"] = error_string
        exceptions = transfers.iloc[exceptions_idx[0] + 2 : details_idx[0]]
        return (
            exceptions[["plate", "row", "column", "error"]]
            .to_numpy(dtype="str")
            .tolist()
        )
    else:
        return []


def read_nexus_repeated_tranfers(path: Path, error_string: str) -> List[str]:
    """
    Import a CSV transfer file and add everything under details to success records.
    Outside of this function, these should be substracted from the error list.
    :param path: path to transfer file (csv)
    :param error_string: verbose description of the error to add to the error list
    :return: List of transfers successful on second try
    """
    transfers = pd.read_csv(
        path,
        header=None,
        names=[
            "Source Plate Name",
            "Source Plate Barcode",
            "Source Plate Type",
            "Source Well",
            "Source Concentration",
            "Source Concentration Units",
            "Destination Plate Name",
            "Destination Plate Barcode",
            "Destination Well",
            "Destination Concentration",
            "Destination Concentration Units",
            "Compound Name",
            "Transfer Volume",
            "Actual Volume",
            "Transfer Status",
            "Current Fluid Height",
            "Current Fluid Volume",
            "% DMSO",
        ],
    )
    if "[DETAILS]" in transfers.values:
        details_idx = transfers.loc[transfers["Source Plate Name"] == "[DETAILS]"].index
        transfers["plate"] = transfers["Destination Plate Barcode"].str.strip(
            "Synthesis"
        )
        transfers["row"] = transfers["Destination Well"].str[0]
        transfers["column"] = transfers["Destination Well"].str[1:]
        transfers["error"] = error_string
        successes = transfers.iloc[details_idx[0] + 2 :]
        return (
            successes[["plate", "row", "column", "error"]]
            .to_numpy(dtype="str")
            .tolist()
        )
    else:
        return []


def read_nexus_survey_errors(
    path: Path, error_string: str, volume_threshold: float
) -> List[str]:
    """
    Import a CSV survey file and add all wells below the threshold volume to error records.
    :param path: path to survey file (csv)
    :param error_string: verbose description of the error to add to the error list
    :param volume_threshold: value below which the volume (in uL) is considered as an error
    :return: List of survey errors
    """
    volumes = pd.read_csv(path, header=7)
    volumes = volumes.drop(columns=["Survey Status"]).dropna(how="any", axis=0)
    volumes["plate"] = volumes["Source Plate Barcode"].str.strip("Synthesis")
    volumes["row"] = volumes["Source Well"].str[0]
    volumes["column"] = volumes["Source Well"].str[1:]
    volumes["error"] = error_string
    exceptions = volumes.loc[
        (volumes["Survey Fluid Volume"] < volume_threshold)
        & (volumes["column"].astype(int).between(3, 22)),
        ["plate", "row", "column", "error"],
    ].astype(object)
    return exceptions.to_numpy(dtype="str").tolist()


def read_mobias_analysis_errors(path: Path, plate_number: int) -> List[str]:
    """
    Import a CSV analysis results file and identify errors in the analysis
    Errors are:
    - too many peaks for target compound or internal standard
    - too strong deviation in the peak area of internal standard
    :param path: path to results file (csv)
    :param plate_number: number of the plate within one experiment
    :return: List of analysis errors
    """

    results = pd.read_csv(path, header=3, encoding="latin-1", skip_blank_lines=False)
    results["plate"] = plate_number
    results["row"] = results["Vial Pos"].apply(lambda x: x.split("-")[1])
    results["column"] = results["Vial Pos"].apply(lambda x: x.split("-")[2])
    results["error_1"] = ""
    results["error_2"] = ""
    results["error_3"] = ""
    results["error_4"] = ""
    internal_standard_number = get_internal_standard_number(path.parent)
    # identify where product A gives multiple peaks. WARNING for 2-4 peaks, ERROR for >4 peaks.
    results.loc[
        (results["SumF1 Cmp"] > 1) & (results["SumF1 Cmp"] <= 4), ["error_1"]
    ] = (
        f"WARNING: multiple peaks for product A ("
        + results["SumF1 Cmp"].astype("str")
        + ")"
    )
    results.loc[(results["SumF1 Cmp"] > 4), ["error_1"]] = (
        f"ERROR: multiple peaks for product A ("
        + results["SumF1 Cmp"].astype(str)
        + ")"
    )
    # identify where IS gives multiple peaks
    results.loc[(results[f"{internal_standard_number} Cmp"] > 1), ["error_2"]] = (
        f"ERROR: multiple peaks for IS ("
        + results[f"{internal_standard_number} Cmp"].astype("str")
        + ")"
    )
    # identify where IS response <50% of plate median
    median_response_area = results[f"{internal_standard_number} Area"].median()
    results.loc[
        results[f"{internal_standard_number} Area"] < median_response_area * 0.5,
        ["error_3"],
    ] = "ERROR: IS response <50% of plate median"

    # identify where IS response >200% of plate median
    results.loc[
        results[f"{internal_standard_number} Area"] > median_response_area * 2.0,
        ["error_4"],
    ] = "ERROR: IS response >200% of plate median"

    # results["error"] = results.apply(aggregate_errors, axis=1)
    results["error"] = (
        results[["error_1", "error_2", "error_3", "error_4"]]
        .agg("; ".join, axis=1)
        .str.strip("; ")
    )
    exceptions = results.loc[
        results["error"] != "", ["plate", "row", "column", "error"]
    ]
    return exceptions.to_numpy(dtype="str").tolist()


def get_nexus_errors(nexus_dir: Path) -> List[str]:
    """
    Iterate through all NEXUS transfer and survey files and collect the errors
    :param nexus_dir: directory under which NEXUS transfer files are stored
    :return: List of transfer errors
    """
    error_list = []
    # task 1: Iterate I_M transfers
    for child in (nexus_dir / "I_M").iterdir():
        if child.is_file() and "transfer" in child.name:
            error_list += read_nexus_transfer_errors(child, "ERROR: I/M transfer error")
    if (nexus_dir / "I_M" / "repeated_transfers").exists():
        for child in (nexus_dir / "I_M" / "repeated_transfers").iterdir():
            if child.is_file() and "transfer" in child.name:
                for successful_transfer in read_nexus_repeated_tranfers(
                    child, "ERROR: I/M transfer error"
                ):
                    error_list.remove(successful_transfer)
    # task 2: Iterate T transfers
    for child in (nexus_dir / "T").iterdir():
        if child.is_file() and "transfer" in child.name:
            error_list += read_nexus_transfer_errors(child, "ERROR: T transfer error")
    if (nexus_dir / "T" / "repeated_transfers").exists():
        for child in (nexus_dir / "T" / "repeated_transfers").iterdir():
            if child.is_file() and "transfer" in child.name:
                for successful_transfer in read_nexus_repeated_tranfers(
                    child, "ERROR: T transfer error"
                ):
                    error_list.remove(successful_transfer)
    # task 3: Iterate dilution transfers
    for child in (nexus_dir / "dilution").iterdir():
        if child.is_file() and "transfer" in child.name:
            error_list += read_nexus_transfer_errors(
                child, "ERROR: Dilution transfer error"
            )
    # task 4: Iterate plate survey before dilution
    for child in (nexus_dir / "dilution").iterdir():
        if child.is_file() and "survey" in child.name:
            error_list += read_nexus_survey_errors(
                child,
                "ERROR: Dilution survey low volume (<2.2 uL)",
                volume_threshold=2.2,
            )
    return error_list


def get_mobias_errors(exp_nr: int) -> List[str]:
    """
    Iterate through all MoBiAS results files and collect the errors
    :param exp_nr: number of the experiment under consideration
    :return: List of transfer errors
    """
    error_list = []
    plate_list_df = pd.read_csv(PLATE_LIST_PATH)
    plate_list = plate_list_df.loc[
        plate_list_df["exp_nr"] == exp_nr,
        ["plate_nr", "lab_journal_number", "results_file_name"],
    ].values.tolist()
    # sanity check: all plate numbers must be unique
    if len(set([i[0] for i in plate_list])) != len([i[0] for i in plate_list]):
        raise RuntimeError("Plate numbers are not unique for this experiment")
    # sanity check: all lab journal numbers must be unique
    if len(set([i[1] for i in plate_list])) != len([i[1] for i in plate_list]):
        raise RuntimeError("Plate numbers are not unique for this experiment")
    for plate_nr, lab_journal_number, results_file_name in plate_list:
        error_list += read_mobias_analysis_errors(
            PLATES_DIR / lab_journal_number / results_file_name, int(plate_nr)
        )
    return error_list


def aggregate_error_list(error_list: list) -> pd.DataFrame:
    """
    Aggregate the error list into a DateFrame that has at most one row per well
    :param error_list: list of errors
    :return: DataFrame with one row per well
    """
    # use intermediate dictionary for the reduction
    error_dict = {}
    for error in error_list:
        try:
            error_dict[f"{error[0]}-{error[1]}-{error[2]}"] += f"; {error[3]}"
        except KeyError:
            error_dict[f"{error[0]}-{error[1]}-{error[2]}"] = f"{error[3]}"
    # transfer dict content to df for convenient CSV writing
    error_df = pd.DataFrame(columns=["well", "errors"], data=error_dict.items())
    error_df["plate"] = error_df["well"].apply(lambda x: x.split("-")[0])
    error_df["row"] = error_df["well"].apply(lambda x: x.split("-")[1])
    error_df["column"] = error_df["well"].apply(lambda x: x.split("-")[2])
    error_df = error_df.sort_values(
        by=["plate", "row", "column"], key=lambda x: pd.to_numeric(x, errors="ignore")
    )
    return error_df.loc[:, ["plate", "row", "column", "errors"]]


def save_errors_to_db(error_df: pd.DataFrame, exp_nr: int):
    """
    Save the errors to the database. This is done by updating the valid column of the experiments table.
    We also set the valid column to NULL error_df does not contain any errors for this record
    :param error_df: DataFrame with errors
    :param exp_nr: Number of the experiment under consideration
    :return: None
    """
    con = sqlite3.connect(DB_PATH)
    cur = con.cursor()
    # Reset the valid column to NULL
    cur.execute(
        "UPDATE experiments SET valid = NULL WHERE exp_nr = ?",
        (exp_nr,),
    )
    # Prepare the data for executemany
    data = [
        (row["errors"], exp_nr, row["plate"], f'{row["row"]}{row["column"]}')
        for _, row in error_df.iterrows()
    ]
    # Update the valid column with the errors using executemany
    cur.executemany(
        "UPDATE experiments SET valid = ? WHERE exp_nr = ? AND plate_nr = ? AND well = ?;",
        data,
    )
    con.commit()
    con.close()
    return


def main(exp_dir: Path, exp_nr: int, skip_nexus=False, skip_confirmation=False):
    """
    Main function to extract the errors from the experiment
    :param exp_dir: Path to the experiment directory
    :param exp_nr: Number of the experiment under consideration
    :param skip_nexus: Whether to skip the NEXUS errors (e.g. for test experiments not performed at NEXUS). Default: False
    :param skip_confirmation: Whether to skip the confirmation step. Default: False
    :return: None
    """
    error_list = []
    error_list += get_manual_error_records(exp_dir / "manual_errorlist.csv")
    if skip_nexus is False:
        error_list += get_nexus_errors(exp_dir / "transfer_files")
    error_list += get_mobias_errors(exp_nr)
    error_df = aggregate_error_list(error_list)
    error_df.to_csv(exp_dir / "extracted_errors.csv", index=False)
    n_errors = len(error_df.loc[error_df["errors"].str.contains("ERROR")])
    n_warnings = len(
        error_df.loc[
            (error_df["errors"].str.contains("WARNING"))
            & ~(error_df["errors"].str.contains("ERROR"))
        ]
    )
    print(
        f"Found {n_errors} wells with errors and further {n_warnings} wells with warnings for experiment #{exp_nr}"
    )
    print(f'A list of errors had been saved to {exp_dir / "extracted_errors.csv"}.')
    if skip_confirmation is False:
        print(
            "Check the list and make any required changes before continuing this program"
        )
        input("Press ENTER to continue...")
    error_df_approved = pd.read_csv(exp_dir / "extracted_errors.csv")
    print("Writing errors to database...")
    save_errors_to_db(error_df_approved, exp_nr)
    print(
        f"Wrote {len(error_df)} error records in database for experiment #{exp_nr}. For all other records from this experiment, the valid column has been reset to NULL."
    )
    return


if __name__ == "__main__":
    # load config
    exp_dir = PLATES_DIR / conf["exp_dir"]

    # parse command line arguments
    parser = argparse.ArgumentParser(
        description="Extract errors from manual error list, NEXUS transfer files and MoBiAS output"
    )
    parser.add_argument(
        "--exp-nr",
        type=int,
        help="Number of the experiment under consideration",
        default=conf["exp_nr"],
    )
    parser.add_argument(
        "--exp-dir", type=Path, help="Path to the experiment directory", default=exp_dir
    )
    parser.add_argument(
        "--skip-nexus",
        action="store_true",
        help="Skip NEXUS transfer file analysis",
        default=conf["errors"]["skip_nexus"],
    )
    parser.add_argument(
        "-y",
        "--yes",
        action="store_true",
        help="Skip user confirmation",
    )
    args = parser.parse_args()

    # run main function
    main(args.exp_dir, args.exp_nr, args.skip_nexus, args.yes)
