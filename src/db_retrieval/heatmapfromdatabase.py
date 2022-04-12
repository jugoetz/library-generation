"""
Retrieve yields from DB, normalize, and plot in heatmaps.


TODO I need a possibility to produce heatmaps for arbitrary data, not just plates.
     Use cases:
      - Compare plate with a set of other syntheses that should be arranged in the same layout as the plate
      - Plot other formats than a plate (e.g. one monomer, all i, all t on the two axes)
     Currently, this is realized in generate_arbitrary_heatmaps.ipynb
"""

import sqlite3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from definitions import PLATES_DIR, DB_PATH
from utils import get_conf

# configuration
# edit config.yaml to change
conf = get_conf()


def read_yields_from_database(db_path, labjournal_nr):
    """Query DB for entries under one lab journal number and arrange the results in a dataframe."""
    # PART 1: Query database to obtain a list of lists
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    query_result = cur.execute('SELECT well, product_A_lcms_ratio, product_B_lcms_ratio, product_C_lcms_ratio, '
                               'product_D_lcms_ratio, product_E_lcms_ratio, product_F_lcms_ratio, product_G_lcms_ratio,'
                               'product_H_lcms_ratio '
                               'FROM main.experiments '
                               'WHERE lab_journal_number = ?;',
                               (labjournal_nr,)
                               ).fetchall()
    con.close()
    # PART 2: Turn the list of lists into a dataframe. Split the 'well' into row and column
    df = pd.DataFrame(data=query_result,
                      columns=['well', 'product_A_lcms_ratio', 'product_B_lcms_ratio', 'product_C_lcms_ratio',
                               'product_D_lcms_ratio', 'product_E_lcms_ratio', 'product_F_lcms_ratio',
                               'product_G_lcms_ratio', 'product_H_lcms_ratio'])
    # a little trick to swap a few wells if necessary bc of minor errors:
    # df['well'] = df['well'].replace({'N14': 'O14',
    #                                        'O14': 'N14',
    #                                        'N16': 'O16',
    #                                        'O16': 'N16',
    #                                        'N18': 'O18',
    #                                        'O18': 'N18',
    #                                        'N20': 'O20',
    #                                        'O20': 'N20',
    #                                        'N22': 'O22',
    #                                        'O22': 'N22',
    #                                        'N24': 'O24',
    #                                        'O24': 'N24',
    #                                        })
    df['row'] = df['well'].str[0]
    df['column'] = df['well'].str[1:]
    return df


def normalize_yields(df, normalization_constant):
    """Divide every value in a product column by a constant"""
    for column in df.columns:
        if column.startswith('product'):
            df[column] = df[column].div(normalization_constant)
    return df


def get_plot(df, product_type, ax=None):
    """Generate a heatmap from dataframe in a matplotlib.pyplot.axis"""
    if ax is None:
        ax = plt.gca()
    # create the heatmap axis
    sns.heatmap(df * 100,
                vmin=0,
                vmax=100,
                annot=True,
                cbar=False,
                cmap=sns.color_palette(
                    'viridis',
                    as_cmap=True,
                ),
                fmt=".1f",
                square=True,
                ax=ax
                )  # this cmap should work for the colorblind

    # style the heatmap
    ax.xaxis.tick_top()  # move column names to top
    ax.tick_params(left=False, top=False)  # remove ticks
    ax.set(xlabel=None, ylabel=None)  # remove axis labels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)  # rotate row labels
    ax.plot()
    ax.set_title(f"Relative yield of product {product_type} by mass trace [%]")
    return ax


def plot_heatmap(df, product_type, save_path, plate_size):
    """Plot an individual plate-sized heatmap from dataframe"""
    if plate_size == 96:
        plt.figure(figsize=(7.5, 5))  # this size works for 96 wells plates
    else:
        plt.figure(figsize=(12, 10))  # this size (maybe) works for 384 well plates

    # obtain heatmap in axis of the current figure
    get_plot(df, product_type)

    # save to file
    fig = plt.gcf()
    fig.savefig(save_path, dpi=100)

    # show the plot
    plt.show()
    return


def plot_heatmap_overview(df, save_path):
    """Plot an overview of all products of one plate in a 4x2 grid from dataframe"""
    # Generate figure and 8 axes inside
    fig, axs = plt.subplots(4, 2, figsize=(21, 35))

    # pass all axes to the get_plot function to have a heatmap generated inside
    for ax, product_type in zip(np.nditer(axs, flags=['refs_ok']), 'ABCDEFGH'):
        ax = get_plot(df[['row', 'column', f'product_{product_type}_lcms_ratio']]
                      .pivot(columns='column', index='row', values=f'product_{product_type}_lcms_ratio')
                      .sort_index(axis=1, key=lambda x: [int(y) for y in x]),
                      product_type,
                      ax.tolist()  # tolist is necessary bc ax is a np.0darray so tolist returns the value
                      )
    fig.tight_layout()

    # save to file and show
    fig.savefig(save_path, dpi=100)
    plt.show()
    return


def plot_experiment_heatmap_from_database(db_path, exp_nr, exp_dir, normalization_constant, plate_size):
    """Main function for plotting yield heatmaps from values stored in 'experiments' database"""
    yields = read_yields_from_database(db_path, exp_nr)
    yields = normalize_yields(yields, normalization_constant)

    # plot an overview with all individual heatmaps in a 4x2 grid
    plot_heatmap_overview(yields, exp_dir / f"heatmap_{exp_nr}_overview.png")

    # plot all heatmaps
    for product_type in 'ABCDEFGH':
        df = yields[['row', 'column', f'product_{product_type}_lcms_ratio']]

        plot_df = df \
            .pivot(columns='column', index='row', values=f'product_{product_type}_lcms_ratio') \
            .sort_index(axis=1, key=lambda x: [int(y) for y in x])  # generate individual df for every plot
        # plot the heatmap
        plot_heatmap(plot_df, product_type, exp_dir / f"heatmap_{exp_nr}_{product_type}.png", plate_size)
    return


if __name__ == '__main__':
    for exp_nr in conf['lab_journal_numbers']:
        print(f'Now plotting {exp_nr}...')
        exp_dir = PLATES_DIR / exp_nr
        plot_experiment_heatmap_from_database(DB_PATH, exp_nr, exp_dir, conf['normalization_constant'],
                                              conf['well_plate_size'])
    print('Finished plotting!')
