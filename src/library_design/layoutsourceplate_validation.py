"""
From an (ORDERED!) list of building blocks, generate source plate layouts for Echo-based synthesis.
This script is only for the 1D validation plates, which have a different layout from the library plates.

Inputs:
    - synthesis_plan.json: JSON of a list of experiments. Every experiment contains 1 list for every of the
      3 building blocks. The type of the list must be List[List[str], List[str], List[str]].
      Typically, this is the output of generatelibraryplan.py

Outputs:
    - (In folders exp{nr.}): source plate layout files source_plate_layout_echo.csv: 1 file per folder,
    the number of folders equals len(input list)
"""
from src.definitions import PLATES_DIR
from labware.plates import Plate384


synthesis_plan = [
    [
        [f"I{i}" for i in range(1, 13)],
        [f"M{i}" for i in range(1, 11)],
        [f"T{i}" for i in range(1, 9)],
    ]
]

for exp_nr, exp in enumerate(synthesis_plan):
    print(f"{exp_nr}  |  {exp}")

for exp_nr, exp in enumerate(synthesis_plan):
    print(f"Experiment Nr. {exp_nr + 1}:")
    print(exp[0])  # I
    print(exp[1])  # M
    print(exp[2])  # T

    # some controls
    assert len(exp[0]) == 12
    assert len(exp[1]) == 10
    assert len(exp[2]) == 8

    # instantiate source plate
    source_plate = Plate384(max_vol=65000, dead_vol=15000)

    # now initiators will fill the row A (A1-A24). Two wells per initiator (left to right, top to bottom)
    # volumes are 63 uL for the first and 50 uL for the last well
    for initiator in exp[0]:
        source_plate.fill_well(source_plate.free(), initiator, 63000)
        source_plate.fill_well(source_plate.free(), initiator, 50000)

    # to be able to use Plate.free() for the Monomers and Terminators as well, we use a little trick and fill the empty
    # rows with placeholders that we will get rid of later
    source_plate.fill_span("B1", "E24", "placeholder", 65000)
    source_plate.fill_span("F21", "F24", "placeholder", 65000)
    source_plate.fill_span("G1", "J24", "placeholder", 65000)

    # now we fill monomers into row F (F1-F20). 2 wells per monomer, again in latin reading order.
    # volumes are 63 uL for the first and 60 uL for the second well
    for monomer in exp[1]:
        source_plate.fill_well(source_plate.free(), monomer, 63000)
        source_plate.fill_well(source_plate.free(), monomer, 60000)

    # now we fill terminators into row K.
    # 3 wells per terminator, again in latin reading order.
    # volumes are 65 uL for the first two and 50 uL for the last well
    for terminator in exp[2]:
        for i in range(2):
            source_plate.fill_well(source_plate.free(), terminator, 65000)
        source_plate.fill_well(source_plate.free(), terminator, 50000)

    # as the last source compound, we fill oxalic acid (X) into the entire bottom row P (more than needed, as fallback)
    source_plate.fill_span("P1", "P24", "X", 65000)

    # finally, we empty the placeholder rows.
    source_plate.empty_span("B1", "E24")
    source_plate.empty_span("F21", "F24")
    source_plate.empty_span("G1", "J24")

    print(source_plate)

    exp_dir = PLATES_DIR / "new" / "exp101"
    exp_dir.mkdir(parents=True, exist_ok=True)
    # print plates to csv files
    source_plate.to_csv(exp_dir / "source_plate_layout.csv", save_volumes=True)
