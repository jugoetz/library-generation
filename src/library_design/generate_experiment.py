"""
!!! EXPERIMENTAL !!!

Design a random SynFerm experiment under explorative setting.

For a 1920-reaction experiment, we randomly draw 16 initiators, 12 monomers, 10 terminators.
To encourage explorative behavior, we query the previously conducted reactions for usage numbers and weight
building blocks so that compounds that have been used less are more likely chosen for the next experiment.

After choosing compounds, we append the experiment to the experiment plan and generate a folder holding relevant
information.

TODO:
    - Query previous experiments
    - Query possible building blocks (i.e. all - missing bc too little material etc.) SIDE QUESRION: What would be needed to pull new building blocks from ChemInventory?
    - Assign weights based on previous usage numbers
    - Weighted sampling
    - Write to experiment plan JSON and CSV
    -
"""
