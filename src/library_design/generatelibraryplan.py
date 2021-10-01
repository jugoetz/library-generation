"""
Generate a randomized Library Synthesis Plan for 50000 compounds.

Input:
    - inventory_compounds.csv: (Pruned) overview of the compounds for the synthesis.
      Taken from cheminventory_cleanup.ipynb.

Output:
    - synthesis_plan.csv: List of experiments to be conducted for the synthesis of 50000 compounds.
                            Every experiment consists of 6 plates (or 1920 wells). For all experiments, a list of
                            the building blocks is given
    - synthesis_plan.json: Serialized list of the syntheses
    - compound_mapping.txt: Space-delimited table of the shorthand name for a compound (e.g. "M1") and the
      longhand equivalent (e.g. Mon002)

"""
import random as rd
import csv
import json
import pandas as pd

from definitions import BUILDING_BLOCKS_DIR, LIB_INFO_DIR, LOG_DIR

"""seed randomizer for reproducible result"""
seed = 42
rd.seed(seed)

"""experiment design parameters"""
initiators_per_run = 16
monomers_per_run = 12
terminators_per_run = 10
plates_per_run = 6
compounds_per_plate = 320
compounds_per_run = plates_per_run * compounds_per_plate
total_runs = 26
total_targets = compounds_per_run * total_runs  # is 49920

"""import the building blocks from (processed) cheminventory export"""
compounds = pd.read_csv(BUILDING_BLOCKS_DIR / 'inventory_compounds.csv')

"""Add a mapping between shorthand and longhand names to the Dataframe"""
initiator_shorts = [f'I{i + 1}' for i in range(len(compounds.loc[compounds['Category'] == 'I']))]
monomer_shorts = [f'M{i + 1}' for i in range(len(compounds.loc[compounds['Category'] == 'M']))]
terminator_shorts = [f'T{i + 1}' for i in range(len(compounds.loc[compounds['Category'] == 'T']))]

compounds['shorts'] = ''
compounds.loc[compounds['Category'] == 'I', ['shorts']] = initiator_shorts
compounds.loc[compounds['Category'] == 'M', ['shorts']] = monomer_shorts
compounds.loc[compounds['Category'] == 'T', ['shorts']] = terminator_shorts

initiators_all = compounds['shorts'].loc[compounds['Category'] == 'I'].tolist()
monomers_all = compounds['shorts'].loc[compounds['Category'] == 'M'].tolist()
terminators_all = compounds['shorts'].loc[compounds['Category'] == 'T'].tolist()

"""Print some information: average usage of each building block"""
initiators_average = total_targets / len(initiators_all)
monomers_average = total_targets / len(monomers_all)
terminators_average = total_targets / len(terminators_all)
print(f'Each initiator will be used in {initiators_average} reactions on average. '
      f'At 1920 reactions per run (6 plates), that equals using it in {initiators_average / monomers_per_run / terminators_per_run} runs')
print(f'Each monomer will be used in {monomers_average} reactions on average. '
      f'At 1920 reactions per run (6 plates), that equals using it in {monomers_average / initiators_per_run / terminators_per_run} runs')
print(f'Each terminator will be used in {terminators_average} reactions on average. '
      f'At 1920 reactions per run (6 plates), that equals using it in {terminators_average / monomers_per_run / initiators_per_run} runs')

"""
Randomize the building block lists.
Ensure the number of building blocks is a multiple of the number that is used per run.
"""
initiators_number = len(initiators_all) // initiators_per_run * initiators_per_run
monomers_number = len(monomers_all) // monomers_per_run * monomers_per_run
terminators_number = len(terminators_all) // terminators_per_run * terminators_per_run
initiators = rd.sample(initiators_all, k=initiators_number)
monomers = rd.sample(monomers_all, k=monomers_number)
terminators = rd.sample(terminators_all, k=terminators_number)
print(f'Using I, M, T: {initiators_number, monomers_number, terminators_number}')
print(f'Using Initiators: {initiators}')
print(f'Using Monomers: {monomers}')
print(f'Using Terminators: {terminators}')

"""split the randomized lists into portions of size <bb>_per_run"""
initiator_sets = [initiators[i: i + initiators_per_run] for i in range(0, len(initiators), initiators_per_run)]
monomers_sets = [monomers[i: i + monomers_per_run] for i in range(0, len(monomers), monomers_per_run)]
terminators_sets = [terminators[i: i + terminators_per_run] for i in range(0, len(terminators), terminators_per_run)]

"""generate all possible combinations of the building block sets"""
product_set = [[i, m, t]
               for i in initiator_sets
               for m in monomers_sets
               for t in terminators_sets
               ]

"""from the total product set, choose the desired number of combinations randomly"""
synthesis_plan = rd.sample(product_set, k=total_runs)
print(f'Synthesis plan:\n{synthesis_plan}')

"""write a log about used and unused building blocks (originating from random selection)"""
with open(LOG_DIR / 'randomization.log', 'w') as file:
    file.write('Used initiators\n')
    for i in initiators:
        file.write(f'{i}\n')
    file.write('\nUnused initiators\n')
    for i in set(initiators_all) - set(initiators):
        file.write(f'{i}\n')
    file.write('\nUsed monomers\n')
    for m in monomers:
        file.write(f'{m}\n')
    file.write('\nUnused monomers\n')
    for m in set(monomers_all) - set(monomers):
        file.write(f'{m}\n')
    file.write('\nUsed terminators\n')
    for t in terminators:
        file.write(f'{t}\n')
    file.write('\nUnused terminators\n')
    for t in set(terminators_all) - set(terminators):
        file.write(f'{t}\n')

"""write the synthesis plan to csv and json"""
with open(LIB_INFO_DIR / 'synthesis_plan.csv', 'w') as file:  # CSV
    writer = csv.writer(file)
    counter = 1
    for run in synthesis_plan:
        writer.writerow([counter, ])
        for building_block in run:
            writer.writerow(building_block)
        counter += 1

with open(LIB_INFO_DIR / 'synthesis_plan.json', 'w') as file:  # json
    json.dump(synthesis_plan, file)

"""write the mapping between short- and longhand names to txt-file"""
compounds.sort_values(by='Category', inplace=True, kind='mergesort')
with open(BUILDING_BLOCKS_DIR / 'compound_mapping.txt', 'w') as file:
    for i, data in compounds[['shorts', 'Compound Name']].iterrows():
        file.write(' '.join([data['shorts'], data['Compound Name']]) + '\n')
