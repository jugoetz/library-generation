"""
From a list of building blocks, generate target plate layouts that provide all possible combinations.

Inputs:
    - synthesis_plan.json: JSON of a list of experiments. Every experiment contains 1 list for every of the
      3 building blocks. The type of the list must be List[List[str], List[str], List[str]].
      Typically this is the output of generatelibraryplan.py

Outputs:
    - Folders exp{nr.} containing plate layout files: The number of folders equals len(input list)
"""
import json
from copy import deepcopy

from src.definitions import LIB_INFO_DIR, PLATES_DIR
from labware.plates import Plate384Echo

with open(LIB_INFO_DIR / "synthesis_plan.json", "r") as file:
    synthesis_plan = json.load(file)

for exp_nr, exp in enumerate(synthesis_plan):
    print(f"{exp_nr}  |  {exp}")

for exp_nr, exp in enumerate(synthesis_plan):
    """
    First, create a "working copy" of x this is necessary because the list comprehension producing synthesis_plan
    gives shallow copies. Using list.pop() would result in that element being removed for all shallow copies and
    would cause the program to fail when it first hits a monomer set it has already seen in a previous exp.
    We prevent this by deep-copying exp before we use it.
    """
    exp = deepcopy(exp)

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

    # with the numbers asserted above, we need 6 plates
    plates = [Plate384Echo() for _ in range(6)]

    for p in plates:
        # assign the initiators
        for row, cmp in zip(p.rows(), exp[0]):
            p.fill_block((row,), p.columns()[2:22], cmp, 1100)
        # assign the monomers
        p.fill_block(p.rows(), p.columns()[2:12], exp[1].pop(0), 1100)
        p.fill_block(p.rows(), p.columns()[12:22], exp[1].pop(0), 1100)
        # assign the terminators
        for col, cmp in zip(p.columns()[2:22], exp[2] * 2):
            p.fill_column(col, cmp, 1100)
        print(p)

    # test if all combinations are unique
    products = []
    for p in plates:
        for row in p.rows():
            for col in p.columns()[2:22]:
                prod = p.compounds(row + col)
                if not prod:
                    print(f"WARNING: Well {row + col} is empty")
                else:
                    products.append(tuple(prod))  # tuple is hashable
    n_prod = len(products)
    if n_prod != 1920:
        print(f"WARNING: {n_prod} product_generator in plate. Expected 1920.")
    if len(set(products)) != len(products):
        print(f"WARNING: Duplicate product_generator detected")

    exp_dir = (
            PLATES_DIR / "new" / f"exp{exp_nr + 1}"
    )  # never write to the main PLATES_DIR to avoid overwriting manual changes
    exp_dir.mkdir(exist_ok=True, parents=True)
    # print plates to csv files
    for i, p in enumerate(plates):
        p.to_csv(exp_dir / f"plate_layout_plate{i + 1}.csv", save_volumes=True)
