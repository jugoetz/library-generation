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
    return wells[cmpd_nr]


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


test_file = PLATES_DIR / "exp13" / "source_plate_layout.csv"
test_output = PLATES_DIR / "exp13" / "ot2_transfers.csv"

target_plate = Plate384(max_vol=65000, dead_vol=15000)
target_plate.from_csv(test_file)

transfers = defaultdict(lambda: ([], []))

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

df[["step", "source_plate", "source_well", "target_wells", "target_volumes"]].to_csv(test_output, index=False)
