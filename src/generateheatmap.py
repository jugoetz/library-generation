#!/usr/bin/env python


"""
This is the 50k project version of generateheatmap, which is tailored to the folder structure in this project.
"""


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from pathlib import Path


DATA_DIR = Path('..', 'data').resolve()
OUTPUT_DIR = DATA_DIR / 'outputs'
INPUT_DIR = DATA_DIR / 'inputs'
EXP_DIR = OUTPUT_DIR / 'target_plates' / 'test_plates_JG213'
DEBUG = True
NORMALIZATION_CONSTANT = 6.0
results_file = EXP_DIR / 'BMII001985_Skript-Results_Variante_B_Lock_322.csv'

well_position_file = EXP_DIR / 'SampleTable_JG213.xls'
norm = 'IS'
IS_compound_number = 8
with open(EXP_DIR / 'notebook_nr.txt') as file:
    LCMS_number = file.read().strip('\n').strip()


def calculate_yield(dict_df, method='IS'):
    # calculate the mass response ratio (i.e. our approximation of yield) for a compounds
    # vs. either internal standard or BPC (i.e. entire area under curve)
    if method == "BPC":
        raise NotImplementedError('Implementation needs to be checked')
        # for k in dict_df.keys():
        #     dict_df[k]["yield"] = dict_df[k]["yield"] / dict_df["BPC"]["yield"]
        # TODO calculate the total yield
        # if type(k) == int:
        #      if k > 1:
        #         df["Total"]["yield"] += df[k]["yield"]
        # print(dict_df)
    elif method == "IS":
        # divide all values by the IS value and normalize (with respect to a hardcoded normalization constant)

        for key, df in dict_df.items():
            df["yield"] = df["yield"] / dict_df["IS"]["yield"] / NORMALIZATION_CONSTANT
            df.fillna(value=0.0, inplace=True)  # TODO this might be a bad idea
    else:
        sys.exit("ERROR: method must be 'BPC' or 'IS', not '{}'".format(method))
    return dict_df


def plot_heatmap(df, measurement_number, method, filename):
    plt.figure(figsize=(7.5, 5))  # this size works for 96 wells plates
    # plt.figure(figsize=(15,10))  # this size (maybe) works for 384 well plates
    axs = sns.heatmap(df * 100,
                      vmin=0,
                      vmax=100,
                      annot=True,
                      cbar=False,
                      cmap=sns.cubehelix_palette(8),
                      fmt=".1f",
                      square=True
                      )  # this cmap should work for the colorblind
    axs.xaxis.tick_top()  # move column names to top
    axs.tick_params(left=False, top=False)  # remove ticks
    axs.set(xlabel=None, ylabel=None)  # remove axis labels
    axs.set_yticklabels(axs.get_yticklabels(), rotation=0)  # rotate row labels
    axs.plot()
    axs.set_title(f"Relative yield of {filename} by mass trace [%]")
    plt.draw()
    plt.savefig(EXP_DIR / f"heatmap_{measurement_number}_{method}_{filename}.png", dpi=200)


# MAIN

if __name__ == '__main__':

    # read data
    results_df = pd.read_csv(results_file, header=3, encoding='latin-1')  # read results.csv file from Mobias
    results_df["Sample ID"] = results_df["Sample ID"].astype("string")  # convert column Sample ID to string
    results_df["Sample ID"] = results_df["Sample ID"].str.split(" ").str[-1]  # extract the JG2xx-001 part
    results_df.columns = results_df.columns.str.strip()  # workaround for inconsistent spaces in Mobias output

    # Form a new df with only the relevant info
    try:
        results_df_clean = results_df[["Sample ID",
                                       "SumF1 Area",
                                       "SumF2 Area",
                                       "SumF3 Area",
                                       "SumF4 Area",
                                       "SumF5 Area",
                                       "SumF6 Area",
                                       "SumF7 Area",
                                       "SumF8 Area",
                                       "BPC Area",
                                       ]]
        print(results_df_clean)
    except KeyError:
        print(f'Available columns: {results_df.columns}')
        raise KeyError

    # Lookup well positions from Sample ID
    wells_df = pd.read_excel(well_position_file)
    wells_df_clean = wells_df[["SampleID.1", "Vial"]].dropna(axis=0)  # TODO the column name might change
    wells_df_clean.rename(columns={'SampleID.1': 'Sample ID'}, inplace=True)
    wells_df_clean = wells_df_clean[
        ~wells_df_clean.iloc[:, 0].str.contains('blank')]  # remove blanks if they are in the input
    wells_df_clean["plate"] = wells_df_clean["Vial"].str.split("-").str[0]
    wells_df_clean["row"] = wells_df_clean["Vial"].str.split("-").str[1]
    wells_df_clean["column"] = wells_df_clean["Vial"].str.split("-").str[2].astype("int")
    wells_df_clean.drop(columns=["Vial"], inplace=True)
    print(wells_df_clean)
    results_merged_df = pd.merge(results_df_clean, wells_df_clean, how="right", on="Sample ID")
    results_merged_df.drop(columns=["Sample ID"], inplace=True)
    print(results_merged_df)

    # Generate individual dataframes for each product of interest
    results = {}
    for i in range(1, len(results_merged_df.columns) - 3):  # -1 each for plate, row, column, BPC, +1 for starting at 1
        results[i] = results_merged_df[['plate', 'row', 'column', f'SumF{i} Area']]
        results[i].rename(columns={f'SumF{i} Area': 'yield'}, inplace=True)
    # Generate dataframe for BPC
    results['BPC'] = results_merged_df[['plate', 'row', 'column', 'BPC Area']].rename(columns={'BPC Area': 'yield'})
    # Generate dataframe for IS
    if norm == 'IS':
        results['IS'] = results[IS_compound_number]
        del results[IS_compound_number]
    print(results)

    # calculate yields from areas
    yields = calculate_yield(results, norm)
    print(yields)

    # plot heatmaps
    for i in yields.keys():
        if i == 'BPC':
            continue  # omit the BPC column when generating plots
        elif i == 'IS':
            continue  # omit the IS column when generating plots
        plot_df = yields[i].pivot(columns='column', index='row',
                                  values='yield')  # generate individual df for every plot
        # plot the heatmap
        plot_heatmap(plot_df, LCMS_number, norm, f'SumF{i}')
        plt.show()  # make sure plots stay open
