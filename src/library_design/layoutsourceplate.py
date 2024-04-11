"""
From an (ORDERED!) list of building blocks, generate source plate layouts for Echo-based synthesis.

Inputs:
    - synthesis_plan.json: JSON of a list of experiments. Every experiment contains 1 list for every of the
      3 building blocks. The type of the list must be List[List[str], List[str], List[str]].
      Typically, this is the output of generatelibraryplan.py

Outputs:
    - (In folders exp{nr.}): source plate layout files source_plate_layout_echo.csv: 1 file per folder,
    the number of folders equals len(input list)
"""
import json

from src.definitions import LIB_INFO_DIR, PLATES_DIR
from labware.plates import Plate384

with open(LIB_INFO_DIR / "synthesis_plan.json", "r") as file:
    synthesis_plan = json.load(file)

for exp_nr, exp in enumerate(synthesis_plan):
    print(f"{exp_nr}  |  {exp}")

for exp_nr, exp in enumerate(synthesis_plan):
    # skip for plates that were already run when the skript was changed
    if exp_nr < 12:  # (starts at 0)
        continue

    print(f"Experiment Nr. {exp_nr + 1}:")
    print(exp[0])  # I
    print(exp[1])  # M
    print(exp[2])  # T

    # some controls
    assert len(exp[0]) == 16  # the 16 KAT will be spread over rows
    assert (
        len(exp[1]) == 12
    )  # the 12 Mon will be spread into 16x10 blocks, giving 2 blocks per plate
    assert len(exp[2]) == 10  # the 10 Ter will be spread over columns

    # instantiate source plate
    source_plate = Plate384(max_vol=65000, dead_vol=15000)

    # now initiators will fill the rows A and B. Three wells per initiator (left to right, top to bottom)
    # volumes are 63 uL for the first two and 35 uL for the last well
    for initiator in exp[0]:
        for i in range(2):
            source_plate.fill_well(source_plate.free(), initiator, 63000)
        source_plate.fill_well(source_plate.free(), initiator, 35000)

    # to be able to use Plate.free() for the Monomers and Terminators as well, we use a little trick and fill the empty
    # rows with placeholders that we will get rid of later
    source_plate.fill_span("C1", "E24", "placeholder", 65000)
    source_plate.fill_span("H1", "J24", "placeholder", 65000)

    # now we fill monomers into rows F and G. 4 wells per monomer, again in latin reading order.
    # volumes are 63 uL for the first three and 25 uL for the last well
    for monomer in exp[1]:
        for i in range(3):
            source_plate.fill_well(source_plate.free(), monomer, 63000)
        source_plate.fill_well(source_plate.free(), monomer, 25000)

    # now we fill terminators into rows K and L and the first 2 wells of M.
    # 5 wells per terminator, again in latin reading order.
    # volumes are 65 uL for the first three and 30 uL for the last well
    for terminator in exp[2]:
        for i in range(4):
            source_plate.fill_well(source_plate.free(), terminator, 65000)
        source_plate.fill_well(source_plate.free(), terminator, 30000)

    # as the last source compound, we fill oxalic acid (X) into the entire bottom row P (more than needed, as fallback)
    source_plate.fill_span("P1", "P24", "X", 65000)

    # finally, we empty the placeholder rows.
    source_plate.empty_span("C1", "E24")
    source_plate.empty_span("H1", "J24")

    print(source_plate)

    exp_dir = PLATES_DIR / "new" / f"exp{exp_nr + 1}"
    # print plates to csv files
    source_plate.to_csv(exp_dir / "source_plate_layout.csv", save_volumes=True)
