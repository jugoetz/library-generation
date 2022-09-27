"""
!!! EXPERIMENTAL, DO NOT USE IN PRODUCTION !!!

Intent:
Design a random SynFerm experiment under explorative setting.

For a 1920-reaction experiment, we randomly draw 16 initiators, 12 monomers, 10 terminators.
To encourage explorative behavior, we query the previously conducted reactions for usage numbers and weight
building blocks so that compounds that have been used less are more likely chosen for the next experiment.

After choosing compounds, we append the experiment to the experiment plan and generate a folder holding relevant
information.

TODO:
    - Query previous experiments
    - Query possible building blocks (i.e. all - missing bc too little material etc.)
       SIDE QUESTION: What would be needed to pull new building blocks from ChemInventory?
    - Assign weights based on previous usage numbers
    - Weighted sampling
    - Write to experiment plan JSON and CSV

Current implementation state:
    ! The current implementation does not actually do any of the above
    - Instead, the current implementation checks how much overlapping reactions we would get
      if we applied the method above
"""
# HOW MUCH OVERLAP WOULD WE GET WITH THIS PARADIGM?

# assume: 4 sets of I, 6 sets of M, 4 sets of T
# I sets containing 16 members, M 12, T 10
# for simplicity, we have 64 I (0-63), 72 M (100-171), 40 T (200-239)
import random
import numpy as np
from copy import copy


def mask_arr(arr, start, end):
    arr = copy(arr)
    arr[:start] = 0
    arr[end:] = 0
    return arr


combinations = []

compounds = list(range(0, 64)) + list(range(100, 172)) + list(range(200, 240))
usage = (
        [3] * 48
        + [0] * 16
        + [4] * 12
        + [2] * 24
        + [1] * 12
        + [0] * 24
        + [6] * 10
        + [2] * 10
        + [1] * 10
        + [0] * 10
)
usage = np.array(usage)
print(usage)

for _ in range(20):
    weights = 1 / (usage + 1) ** 3
    idx_i = random.choices(np.arange(len(compounds)), mask_arr(weights, 0, 64), k=16)
    idx_m = random.choices(np.arange(len(compounds)), mask_arr(weights, 64, 136), k=12)
    idx_t = random.choices(np.arange(len(compounds)), mask_arr(weights, 136, 176), k=10)

    usage[idx_i] += 1
    usage[idx_m] += 1
    usage[idx_t] += 1

    combinations += [
        f"{compounds[i]}_{compounds[m]}_{compounds[t]}"
        for i in idx_i
        for m in idx_m
        for t in idx_t
    ]
print(usage)
print(len(combinations))
print(len(set(combinations)))
