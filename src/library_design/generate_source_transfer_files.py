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
from definitions import PLATES_DIR

EXPERIMENT_FOLDER = "exp13"


def compound_to_well(cmpd):
    """
    Translate by making this substitution i-th compound -> i-th well in a 96 well plate.

    Example:
        compound I15 -> well B3
    """
    if cmpd == "X":  # special case
        return "A1"
    cmpd_nr = int(cmpd[1:])
    wells = [f"{row}{column}" for row in "ABCDEFGH" for column in range(1, 13)]
    return wells[cmpd_nr - 1]  # -1 b/c compound nrs start at 1, but wells is 0-based


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


source_plate_layout_file = PLATES_DIR / EXPERIMENT_FOLDER / "source_plate_layout.csv"
transfer_file = PLATES_DIR / EXPERIMENT_FOLDER / "ot2_transfers.csv"

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

df = pd.DataFrame.from_dict(transfers, orient="index").reset_index().rename(
    columns={"index": "compound", 0: "target_wells", 1: "target_volumes"})

df["source_well"] = df["compound"].apply(compound_to_well)
df["source_plate"] = df["compound"].apply(compound_to_source_plate)
df["step"] = df["source_plate"].apply(lambda x: 2 if x == "terminator" else 1)

# set up transfers for LCMS analysis plate
# logic: for all building blocks -> 1 uL into next well of lcms plate,
# where next well goes in the order A1->D1->A8->D8, and counter is reset for terminators
# practically, we only need to copy the dataframe for I,M,T transfers and adjust target wells and volumes
df_lcms = df.loc[df["source_plate"] != "oxalic acid"].copy()
df_lcms["source_plate"] += "_lcms"
# chained sort to get I1->I99->M1->M99->T1->T99 order
df_lcms = df_lcms.sort_values(
    by="compound", kind="stable", key=lambda x: [int(y[1:]) for y in x]).sort_values(
    by="compound", kind="stable", key=lambda x: [y[0] for y in x])

df_lcms.loc[df_lcms["source_plate"] == "initiator_lcms", "target_wells"] = [f"{r}{c}" for c in range(1, 5) for r in
                                                                            "ABCD"]
df_lcms.loc[df_lcms["source_plate"] == "monomer_lcms", "target_wells"] = [f"{r}{c}" for c in range(5, 8) for r in
                                                                          "ABCD"]
df_lcms.loc[df_lcms["source_plate"] == "terminator_lcms", "target_wells"] = [f"{r}{c}" for c in range(1, 4) for r in
                                                                             "ABCD"][:-2]
df_lcms["target_volumes"] = [1000, ] * len(df_lcms)

df = pd.concat((df, df_lcms), ignore_index=True)

df[["step", "source_plate", "source_well", "target_wells", "target_volumes"]].to_csv(transfer_file, index=False)
