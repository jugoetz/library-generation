"""
Generates the transfer files for OT2 to prepare source plates that can be used at NEXUS.

This script reads the source_plate_layout for an experiment and produces a CSV file holding the transfers

Each row in the csv contains
step | source_plate | source_well | target_wells | target_volumes

where
- step is 1 for initiator/monomer/oxalic acid and 2 for terminator
- source plate is any of {"initiator", "monomer", "terminator", "oxalic acid"}
- source well is the well of the compound in the source plate (e.g. "A1")
- target_wells is a list of wells in the target plate where the compound should be dispensed (e.g. ["A1", "B1"])
- target_volumes is a list of the same length as target_wells indicating the volume to be dispensed into each
    target well
"""
from collections import defaultdict

import pandas as pd

from labware.plates import Plate384
from src.definitions import PLATES_DIR
from src.util.utils import get_conf

# configuration
# edit config.yaml to change
conf = get_conf()


# option source_plate_layout:
# - "outer_wells": source plate has all compounds in A/H rows (for better visibility)
#   the order is such it increases A1->A12->H1->H12, but we leave wells free to avoid cross contamination.
#   Precisely, the 16 initiators occupy A1-A9, A11, and all odd wells in row H,
#   the 12 monomers fill all odd wells in rows A and H
#   the 10 terminators fill all odd wells in row A and all odd wells in H1-H8
# - "canonical_order": compound I10 is in the 10th well (i.e. A10) and so on


def compound_to_well(cmpd):
    """
    Translate by making this substitution: i-th compound -> i-th well in a 96 well plate.

    Example:
        compound I15 -> well B3
    """
    if cmpd == "X":  # special case
        return "A1"
    cmpd_nr = int(cmpd[1:])
    wells = [f"{row}{column}" for row in "ABCDEFGH" for column in range(1, 13)]
    return wells[cmpd_nr - 1]  # -1 b/c compound nrs start at 1, but wells is 0-based


def compound_to_outer_well(cmpds, cmpd_type):
    cmpd_list = sorted(cmpds.values.tolist(), key=lambda x: (x[0], int(x[1:])))
    if cmpd_type == "I":
        wells = [
            "A1",
            "A2",
            "A3",
            "A4",
            "A5",
            "A6",
            "A7",
            "A8",
            "A9",
            "A11",
            "H1",
            "H3",
            "H5",
            "H7",
            "H9",
            "H11",
        ]
    elif cmpd_type == "M":
        wells = [
            "A1",
            "A3",
            "A5",
            "A7",
            "A9",
            "A11",
            "H1",
            "H3",
            "H5",
            "H7",
            "H9",
            "H11",
        ]
    elif cmpd_type == "T":
        wells = ["A1", "A3", "A5", "A7", "A9", "A11", "H1", "H3", "H5", "H7"]
    else:
        raise ValueError("Invalid cmpd_type. Options: 'I', 'M', 'T'")
    well_map = dict(zip(cmpd_list, wells))
    return cmpds.map(well_map)


def compound_to_source_plate(cmpd):
    """
    Translate compound to source plate

    Example:
        compound I15 -> initiator
    """
    if cmpd == "X":
        return "oxalic acid"
    elif cmpd.startswith("I"):
        return "initiator"
    elif cmpd.startswith("M"):
        return "monomer"
    elif cmpd.startswith("T"):
        return "terminator"


source_plate_layout_file = PLATES_DIR / conf["exp_dir"] / "source_plate_layout.csv"
transfer_file = PLATES_DIR / conf["exp_dir"] / "ot2_transfers.csv"

print(f"Reading source plate layout from {source_plate_layout_file}...")

target_plate = Plate384(max_vol=65000, dead_vol=15000)
target_plate.from_csv(source_plate_layout_file)

transfers = defaultdict(lambda: ([], []))

# set up transfers from target plate layout
for well in target_plate.wells():
    well_content = target_plate.well(well)
    if well_content[1] == 0:
        continue
    transfers[well_content[0][0]][0].append(well)
    transfers[well_content[0][0]][1].append(well_content[1])

df = (
    pd.DataFrame.from_dict(transfers, orient="index")
    .reset_index()
    .rename(columns={"index": "compound", 0: "target_wells", 1: "target_volumes"})
)
df["source_plate"] = df["compound"].apply(compound_to_source_plate)
if conf["ot2_transfers"]["source_plate_layout"] == "canonical_order":
    df["source_well"] = df["compound"].apply(compound_to_well)
elif conf["ot2_transfers"]["source_plate_layout"] == "outer_wells":
    df[
        "source_well"
    ] = "A1"  # by initializing as A1, we don't have to take care of oxalic acid wells
    df.loc[df["source_plate"] == "initiator", "source_well"] = compound_to_outer_well(
        df.loc[df["source_plate"] == "initiator", "compound"], "I"
    )
    df.loc[df["source_plate"] == "monomer", "source_well"] = compound_to_outer_well(
        df.loc[df["source_plate"] == "monomer", "compound"], "M"
    )
    df.loc[df["source_plate"] == "terminator", "source_well"] = compound_to_outer_well(
        df.loc[df["source_plate"] == "terminator", "compound"], "T"
    )

else:
    raise ValueError("Invalid option for 'source_plate_layout'.")

df["step"] = df["source_plate"].apply(lambda x: 2 if x == "terminator" else 1)

# set up transfers for LCMS analysis plate
# logic: for all building blocks -> 1 uL into next well of lcms plate,
# where next well goes in the order A1->D1->A8->D8, and counter is reset for terminators
# practically, we only need to copy the dataframe for I,M,T transfers and adjust target wells and volumes
df_lcms = df.loc[df["source_plate"] != "oxalic acid"].copy()
df_lcms["source_plate"] += "_lcms"
# chained sort to get I1->I99->M1->M99->T1->T99 order
df_lcms = df_lcms.sort_values(
    by="compound", kind="stable", key=lambda x: [int(y[1:]) for y in x]
).sort_values(by="compound", kind="stable", key=lambda x: [y[0] for y in x])

df_lcms.loc[df_lcms["source_plate"] == "initiator_lcms", "target_wells"] = [
    f"{r}{c}" for c in range(1, 5) for r in "ABCD"
]
df_lcms.loc[df_lcms["source_plate"] == "monomer_lcms", "target_wells"] = [
    f"{r}{c}" for c in range(5, 8) for r in "ABCD"
]
df_lcms.loc[df_lcms["source_plate"] == "terminator_lcms", "target_wells"] = [
    f"{r}{c}" for c in range(1, 4) for r in "ABCD"
][:-2]
df_lcms["target_volumes"] = [
    1000,
] * len(df_lcms)

df = pd.concat((df, df_lcms), ignore_index=True)

print(f"Writing transfers to {transfer_file}...")

df[["step", "source_plate", "source_well", "target_wells", "target_volumes"]].to_csv(
    transfer_file, index=False
)

print("Done.")
