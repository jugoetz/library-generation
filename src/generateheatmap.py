#!/usr/bin/env python


"""
This is the 50k project version of generateheatmap, which is tailored to the folder structure in this project.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from config import *
import json
import sqlite3
import numpy as np

"""GLOBALS"""
DEBUG = True
NORMALIZATION_CONSTANT = 6.0
results_file = EXP_DIR / 'BMII001987_Skript-Results_protecting-groups.csv'
well_position_file = EXP_DIR / 'BMIIyyyyyy-SampleTable_JG216-big.xls'
norm = 'IS'
stop_after_db_write = True
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
            df["yield"] = df["yield"].astype('float64') / dict_df["IS"]["yield"]
            # TODO have something here to save values to DB before normalization. Also, assign which de-PG product corresponds to which main product
            df.fillna(value=0.0, inplace=True)  # TODO this might be a bad idea bc it hides errors
    else:
        sys.exit("ERROR: method must be 'BPC' or 'IS', not '{}'".format(method))
    return dict_df


def save_to_db(dict_df):
    con = sqlite3.connect(DB_PATH)
    cur = con.cursor()

    def add(df_list):
        """
        Function to add a list of dataframes where each df has 4 columns (plate, row, column, yield) and where the
        first three are used for identification while the 4th is added. The return value is thus a 4-column dataframe.
        The first dataframe must have all the identifiers (first 3 columns) that will be present.
        """
        df = df_list[0]
        for other in df_list[1:]:
            df = df.merge(other, how='left', on=['plate', 'row', 'column'])
        df['cumulated_yield'] = df[[i for i in df.columns if 'yield' in i]].sum(axis=1)
        df = df[['plate', 'row', 'column', 'cumulated_yield']]
        df.rename(columns={'cumulated_yield': 'yield'}, inplace=True)
        return df

    # load the saved connections between product type (A) and SumFx field
    with open(EXP_DIR / 'compound_alternative_mass_dict.json', 'r') as jfile:
        side_product_dict = json.load(jfile)
    # reformat dict for easier use
    side_product_associations = {t: [v for k, v in side_product_dict.items() if k.startswith(t)] for t in 'ABCDEFGH'}
    for k, v in side_product_associations.items():
        product_yields = add([dict_df[int(val.strip('SumF'))] for val in v])
        print(product_yields)
        for i, data in product_yields.iterrows():
            cur.execute(
                f'UPDATE main.experiments SET product_{k}_lcms_ratio = ? WHERE lab_journal_number = ? AND well = ?;',
                (data['yield'], EXP_NR, f'{data["row"]}{data["column"]}'))
        con.commit()
    con.close()
    return


def normalize_yields(dict_df, normalization_constant):
    for df in dict_df.values():
        df['yield'] = df['yield'].div(normalization_constant)
    return dict_df


def plot_heatmap(df, measurement_number, method, filename):
    if PLATE_SIZE == 96:
        plt.figure(figsize=(7.5, 5))  # this size works for 96 wells plates
    else:
        plt.figure(figsize=(12, 10))  # this size (maybe) works for 384 well plates

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
    plt.savefig(EXP_DIR / f"heatmap_{measurement_number}_{method}_{filename}.png", dpi=100)


# MAIN

if __name__ == '__main__':

    # read data
    results_df = pd.read_csv(results_file, header=3, encoding='latin-1')  # read results.csv file from Mobias
    results_df["Sample ID"] = results_df["Sample ID"].astype("string")  # convert column Sample ID to string
    results_df["Sample ID"] = results_df["Sample ID"].str.split(" ").str[-1]  # extract the JG2xx-001 part
    results_df.columns = results_df.columns.str.strip()  # workaround for inconsistent spaces in Mobias output
    results_df.replace(to_replace='-', value=0,
                       inplace=True)  # replace '-' value that is present when no compound to search for was given
    # mind that columns _Area in results_df can have object dtype if '-' was present. Use .astype('int64') downstream when using this for calculations
    # Form a new df with only the relevant info
    cols = results_df.columns
    cols = [s for s in cols if "Area" in s and "UV" not in s and results_df[s].count() > 0]
    cols.insert(0, 'Vial Pos')
    cols.insert(0, 'Sample ID')
    try:
        results_df_clean = results_df[cols]
        print(results_df_clean)
    except KeyError:
        print(f'Available columns: {results_df.columns}')
        raise KeyError

    use_sample_table = False

    if use_sample_table is True:
        # Lookup well positions from Sample ID
        wells_df = pd.read_excel(well_position_file)
        wells_df_clean = wells_df[['Sample-Ident', 'Vial']].dropna(axis=0)  # TODO the column name might change
        wells_df_clean.rename(columns={'Sample-Ident': 'Sample ID'}, inplace=True)
        wells_df_clean = wells_df_clean[
            ~wells_df_clean.iloc[:, 0].str.contains('blank', case=False)]  # remove blanks if they are in the input
        wells_df_clean['plate'] = wells_df_clean['Vial'].str.split('-').str[0]
        wells_df_clean['row'] = wells_df_clean['Vial'].str.split('-').str[1]
        wells_df_clean['column'] = wells_df_clean['Vial'].str.split('-').str[2].astype('int')
        wells_df_clean.drop(columns=['Vial'], inplace=True)
        print(wells_df_clean)
        results_merged_df = pd.merge(results_df_clean, wells_df_clean, how='right', on='Sample ID')
        results_merged_df.drop(columns=['Sample ID'], inplace=True)
        print(results_merged_df)
    else:
        # get well positions from "Vial Pos" (introduced around JG215)
        results_df_clean['plate'] = results_df_clean['Vial Pos'].str.split('-').str[0]
        results_df_clean['row'] = results_df_clean['Vial Pos'].str.split('-').str[1]
        results_df_clean['column'] = results_df_clean['Vial Pos'].str.split('-').str[2]
        results_merged_df = results_df_clean
        results_merged_df.drop(columns=['Sample ID', 'Vial Pos'], inplace=True)

    # Generate individual dataframes for each product of interest
    results = {}
    for i in range(1, len(results_merged_df.columns) - 3):  # -1 each for plate, row, column, BPC, +1 for starting at 1
        results[i] = results_merged_df[['plate', 'row', 'column', f'SumF{i} Area']]
        results[i].rename(columns={f'SumF{i} Area': 'yield'}, inplace=True)
    # Generate dataframe for BPC
    results['BPC'] = results_merged_df[['plate', 'row', 'column', 'BPC Area']].rename(columns={'BPC Area': 'yield'})
    # Generate dataframe for IS
    if norm == 'IS':
        with open(EXP_DIR / 'compound_alternative_mass_dict.json', 'r') as jfile:
            side_product_dict = json.load(jfile)
        IS_compound_number = int(side_product_dict['IS_formula'].strip('SumF'))
        results['IS'] = results[IS_compound_number]
        del results[IS_compound_number]
    print(results)

    # calculate yields from areas
    yields = calculate_yield(results, norm)
    save_to_db(yields)
    if stop_after_db_write is True:
        exit(0)
    yields = normalize_yields(yields, NORMALIZATION_CONSTANT)
    print(yields)

    # plot heatmaps
    for i, df in yields.items():
        if i == 'BPC':
            continue  # omit the BPC column when generating plots
        elif i == 'IS':
            continue  # omit the IS column when generating plots
        plot_df = df \
            .pivot(columns='column', index='row', values='yield') \
            .sort_index(axis=1, key=lambda x: [int(y) for y in x])  # generate individual df for every plot
        # plot the heatmap
        plot_heatmap(plot_df, LCMS_number, norm, f'SumF{i}')
        plt.show()  # make sure plots stay open
