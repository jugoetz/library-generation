"""
Generate the weigh-in sheet for an experiment.
The experiment can be defined by either an exp_nr or a lab_journal_number (both as in the experiments DB table).
    -> To use this with (single plate) non-canonical reactions, use a lab_journal_number
The directory connected to the experiment must be specified (exp_dir). This is where output will be saved.
Volumes can be adjusted in configuration.
"""

from pathlib import Path

import xlsxwriter

from src.definitions import PLATES_DIR
from src.util.db_utils import SynFermDatabaseConnection
from src.util.utils import get_conf

# configuration
# edit config.yaml to change
conf = get_conf()


def main():
    exp_dir = conf["exp_dir"]
    print(f"Generating weigh-in for {exp_dir}...")
    mycon = SynFermDatabaseConnection()

    # get the list of compounds of interest
    if "exp_nr" in conf and "lab_journal_number" in conf["weight"]:
        raise ValueError(
            'Exactly one of "exp_nr" and "weight.lab_journal_number" must be None.'
        )
    elif "exp_nr" in conf:
        compounds = mycon.get_starting_materials_for_experiment(exp_nr=conf["exp_nr"])
    elif "lab_journal_number" in conf["weight"]:
        compounds = mycon.get_starting_materials_for_experiment(
            lab_journal_number=conf["lab_journal_nr"]
        )
    else:
        raise ValueError(
            'Only (and exactly) one of "exp_nr" and "lab_journal_nr" may be None.'
        )

    # for all compounds (shorts), assign long names and needed solution volumes
    compound_volumes = {}
    for i in compounds[0]:  # initiators
        compound_volumes[i] = [
            mycon.get_long_name(i),
            conf["weight"]["initiator_volume"],
        ]
    for m in compounds[1]:  # monomers
        compound_volumes[m] = [mycon.get_long_name(m), conf["weight"]["monomer_volume"]]
    for t in compounds[2]:  # terminators
        compound_volumes[t] = [
            mycon.get_long_name(t),
            conf["weight"]["terminator_volume"],
        ]

    # add weights needed for the minimum solution volume
    for k, v in compound_volumes.items():
        mol_wt = mycon.get_molecular_weight(k)
        # m = M (g/mol) * c (mol/L) * V (uL) / 1000 ug/mg
        mass = mol_wt * 0.05 * v[1] / 1000.0
        v.append(mass)
        print(f"{k}: {v}")

    # flatten dict to list for convenience with the output
    cmp_list = [[k] + v for k, v in compound_volumes.items()]

    # check if output file exists and prompt user if so
    if Path(PLATES_DIR / exp_dir / "weigh-in.xlsx").exists():
        go_on = input(
            f'{PLATES_DIR / exp_dir / "weigh-in.xlsx"} already exists? Overwrite? [[y]/n]:\n'
        )
        if go_on == "y":
            pass
        else:
            exit(1)

    # write the excel sheet
    workbook = xlsxwriter.Workbook(PLATES_DIR / exp_dir / "weigh-in.xlsx")
    worksheet = workbook.add_worksheet()
    decimal_format_one_place = workbook.add_format()
    decimal_format_one_place.set_num_format("#.0")
    decimal_format_int = workbook.add_format()
    decimal_format_int.set_num_format(1)

    worksheet.write(0, 0, "Short")
    worksheet.write(0, 1, "Long")
    worksheet.write(0, 2, "Barcode")
    worksheet.write(0, 3, "theor. m [mg]")
    worksheet.write(0, 4, "actual m [mg]")
    worksheet.write(0, 5, "V [uL]")
    worksheet.write(0, 6, "actual V [uL]")
    worksheet.write(0, 7, "V DMSO [uL]")
    worksheet.write(0, 8, "V oxalic [uL]")
    worksheet.write(0, 9, "comment")
    col = 0
    for i, cmp in enumerate(cmp_list):
        row = i + 1  # +1 for header
        a = cmp[0]
        b = cmp[1]
        f = cmp[2]
        d = round(cmp[3], 2)

        g = f"=E{row + 1}/D{row + 1}*F{row + 1}"  # +1 for excel convention to start at 1
        h = f"=G{row + 1}*0.9"
        i = f"=G{row + 1}*0.1"

        worksheet.write(row, col, a)
        worksheet.write(row, col + 1, b)
        worksheet.write(row, col + 3, d, decimal_format_one_place)
        worksheet.write(row, col + 4, None, decimal_format_one_place)
        worksheet.write(row, col + 5, f, decimal_format_int)
        worksheet.write(row, col + 6, g, decimal_format_int)
        worksheet.write(row, col + 7, h, decimal_format_int)
        worksheet.write(row, col + 8, i, decimal_format_int)

    workbook.close()


if __name__ == "__main__":
    main()
