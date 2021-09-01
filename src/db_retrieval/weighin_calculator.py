import pickle as pkl
import json
import xlsxwriter
from config import *

"""
Right now has many hardcoded assumptions that make it useful only for the 50 k project.
Note that this works on exp{x}, not JGxxx, so this exp_nr needs to be an experiment nr not a lab journal number
"""

# control variable TODO import from a config file
exp_nr = 6


def main(exp_nr):
    exp_dir = PLATES_DIR / f'exp{exp_nr}'

    # import the short-long name relation TODO this can be fetched from DB
    mapping = {}
    with open(BB_DIR / 'compound_mapping.txt', 'r') as file:
        for l in file.readlines():
            mapping[l.split()[0]] = l.split()[1]

    # import properties for weigh-in mass TODO this can be fetched from DB
    with open(LIB_INFO_DIR / 'library_constituents_dataframe.pkl', 'rb') as file:
        df = pkl.load(file)

    # import which compounds are of interest TODO this can be fetched from DB
    with open(LIB_INFO_DIR / 'synthesis_plan.json', 'r') as file:
        synthesis_plan = json.load(file)

    print(f'Generating weigh-in for experiment {exp_nr}...')
    compounds = synthesis_plan[exp_nr - 1]

    # overwrite inputs
    # compounds = [['I15', 'I56'], [f'M{i+1}' for i in range(74)], ['T19', 'T31']]

    i_volume = 2 * 65 + 35  # in uL. this is wells needed * max well volume (last well does not need to be full)
    m_volume = 3 * 65 + 25
    t_volume = 4 * 65 + 30

    # overwrite volumes

    # i_volume = 2 * 65  # in uL. this is wells needed * max well volume
    # m_volume = 1 * 65
    # t_volume = 2 * 65

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
    cmp_list = [[k] + v for k, v in compound_volumes.items()]

    # output
    if Path(exp_dir / 'weigh-in.xlsx').exists():
        go_on = input(f'{exp_dir / "weigh-in.xlsx"} already exists? Overwrite? [[y]/n]:\n')
        if go_on == 'y':
            pass
        else:
            exit(1)

    workbook = xlsxwriter.Workbook(exp_dir / 'weigh-in.xlsx')
    worksheet = workbook.add_worksheet()
    decimal_format_one_place = workbook.add_format()
    decimal_format_one_place.set_num_format('#.0')
    decimal_format_int = workbook.add_format()
    decimal_format_int.set_num_format(1)

    worksheet.write(0, 0, 'Short')
    worksheet.write(0, 1, 'Long')
    worksheet.write(0, 2, 'Barcode')
    worksheet.write(0, 3, 'theor. m [mg]')
    worksheet.write(0, 4, 'actual m [mg]')
    worksheet.write(0, 5, 'V [uL]')
    worksheet.write(0, 6, 'actual V [uL]')
    worksheet.write(0, 7, 'V DMSO [uL]')
    worksheet.write(0, 8, 'V oxalic [uL]')
    worksheet.write(0, 9, 'comment')
    col = 0
    for i, cmp in enumerate(cmp_list):
        row = i + 1  # +1 for header
        a = cmp[0]
        b = cmp[1]
        f = cmp[2]
        d = round(cmp[3], 2)

        g = f'=E{row + 1}/D{row + 1}*F{row + 1}'  # +1 for excel convention to start at 1
        h = f'=G{row + 1}*0.9'
        i = f'=G{row + 1}*0.1'

        worksheet.write(row, col, a)
        worksheet.write(row, col + 1, b)
        worksheet.write(row, col + 3, d, decimal_format_one_place)
        worksheet.write(row, col + 4, None, decimal_format_one_place)
        worksheet.write(row, col + 5, f, decimal_format_int)
        worksheet.write(row, col + 6, g, decimal_format_int)
        worksheet.write(row, col + 7, h, decimal_format_int)
        worksheet.write(row, col + 8, i, decimal_format_int)

    workbook.close()


if __name__ == '__main__':
    main(exp_nr=exp_nr)
