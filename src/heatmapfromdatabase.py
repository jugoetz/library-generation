#!/usr/bin/env python


"""
This is the latter part of the former "generateheatmap", here we read the values from the database. Yields in DB are not normalized
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from config import *
import json
import sqlite3
import numpy as np

NORMALIZATION_CONSTANT = 6.0


def read_yields_from_database(db_path, labjournal_nr):
    # PART 1: Query database to obtain a list of lists
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    query_result = cur.execute(
        'SELECT well, product_A_lcms_ratio, product_B_lcms_ratio, product_C_lcms_ratio, product_D_lcms_ratio, product_E_lcms_ratio, product_F_lcms_ratio, product_G_lcms_ratio, product_H_lcms_ratio FROM main.experiments WHERE lab_journal_number = ?;',
        (labjournal_nr,)).fetchall()
    con.close()
    # PART 2: Turn the list of lists into a dataframe. Split the 'well' into row and column
    df = pd.DataFrame(data=query_result,
                      columns=['well', 'product_A_lcms_ratio', 'product_B_lcms_ratio', 'product_C_lcms_ratio',
                               'product_D_lcms_ratio', 'product_E_lcms_ratio', 'product_F_lcms_ratio',
                               'product_G_lcms_ratio', 'product_H_lcms_ratio'])
    df['row'] = df['well'].str[0]
    df['column'] = df['well'].str[1:]
    return df


def normalize_yields(df, normalization_constant):
    for column in df.columns:
        if column.startswith('product'):
            df[column] = df[column].div(normalization_constant)
    return df


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
    yields = read_yields_from_database(DB_PATH, EXP_NR)
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
        plot_heatmap(plot_df, EXP_NR, 'IS', f'SumF{i}')
        plt.show()  # make sure plots stay open
