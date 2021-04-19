from labware.plates import Plate96, Plate384
from pathlib import Path
import pickle as pkl
import re
import os
import xlsxwriter

"""
Right now has many hardcoded assumptions that make it useful only for the 50 k project.
"""


"""GLOBALS"""
# directories and files
DATA_DIR = Path('..', 'data').resolve()
OUTPUT_DIR = DATA_DIR / 'outputs'
INPUT_DIR = DATA_DIR / 'inputs'
EXP_DIR = OUTPUT_DIR / 'target_plates' / 'exp1'
# PLATE_REGEX = re.compile('test_JG([0-9]+).csv')
PLATE_REGEX = re.compile('plate_layout_plate([0-9]+).csv')
COMPOUND_MAPPING = OUTPUT_DIR / 'compound_mapping.txt'


with open(OUTPUT_DIR / 'library_constituents_dataframe.pkl', 'rb') as file:
    df = pkl.load(file)

with open(OUTPUT_DIR / 'synthesis_plan.pkl', 'rb') as file:
    synthesis_plan = pkl.load(file)

mapping = {}
with open(OUTPUT_DIR / 'compound_mapping.txt', 'r') as file:
    for l in file.readlines():
        mapping[l.split()[0]] = l.split()[1]

exp_nr = int(EXP_DIR.name.strip('exp')) - 1
print(f'Generating weigh-in for experiment {exp_nr}...')
compounds = synthesis_plan[exp_nr]

i_volume = 3 * 65  # in uL. this is wells needed * max well volume
m_volume = 4 * 65
t_volume = 5 * 65

compound_volumes = {}
for i in compounds[0]:  # initiators
    long = mapping[i]
    compound_volumes[i] = [long, i_volume]
for m in compounds[1]:  # monomers
    long = mapping[m]
    compound_volumes[m] = [long, m_volume]
for t in compounds[2]:  # terminators
    long = mapping[t]
    compound_volumes[t] = [long, t_volume]

# add minimum weights
for k, v in compound_volumes.items():
    print(f'{k}: {v}')
    mass_per_100uL = df.loc[df['Compound Name'] == v[0], 'weigh-in [mg] / 100 µL'].values[0]
    mass = mass_per_100uL * v[1] / 100
    v.append(mass)

# flatten dict to list for convenience with the output
cmp_list = [[k]+v for k, v in compound_volumes.items()]
# output
workbook = xlsxwriter.Workbook(EXP_DIR / 'weigh-in.xlsx')
worksheet = workbook.add_worksheet()
worksheet.write(0, 0, 'Short')
worksheet.write(0, 1, 'Long')
worksheet.write(0, 2, 'V [uL]')
worksheet.write(0, 3, 'theor. m [mg]')
worksheet.write(0, 4, 'actual m [mg]')
worksheet.write(0, 5, 'actual V [uL]')
col = 0
for i, cmp in enumerate(cmp_list):
    row = i + 1  # +1 for header
    a = cmp[0]
    b = cmp[1]
    c = cmp[2]
    d = round(cmp[3], 2)

    f = f'=E{row + 1}/D{row + 1}*C{row + 1}'  # +1 for excel convention to start at 1

    worksheet.write(row, col, a)
    worksheet.write(row, col + 1, b)
    worksheet.write(row, col + 2, c)
    worksheet.write(row, col + 3, d)
    worksheet.write(row, col + 5, f)


workbook.close()




