"""
Go from raw lcms analysis results (i.e. peak areas) to a yield metric (i.e. peak area / IS area for every compound).

We start with a DB table 'lcms':
id | synthesis id  | [lcms compounds] | [lcms areas]

We go to a DB table 'lcms yields':
id (foreign key from lcms table) | SumF.. | type (e.g. A) | SMILES (or better: id of the product in a different table)
# TODO what is the best way here? mulitple entries per experiment with
# TODO looks like I need a products table
"""

import json
import sqlite3
import re
import pandas as pd
from config import *


def split_products_into_dataframes(df, exp_dir):
    """Generate individual dataframe for each product (and for special cases like IS and BPC). Return a dict of those"""
    df_dict = {}

    # find the product columns in df by regex
    regex = re.compile('SumF([0-9]+) Area')
    product_columns = filter(regex.match, df.columns)

    # for every product column, add a separate df to dict
    for c in product_columns:
        number = int(regex.match(c).group(1))
        df_dict[number] = df[['plate', 'row', 'column', c]]
        df_dict[number].rename(columns={c: 'yield'}, inplace=True)

    # Add df for BPC to dict as well
    df_dict['BPC'] = df[['plate', 'row', 'column', 'BPC Area']].rename(columns={'BPC Area': 'yield'})

    # Generate dataframe for IS TODO bad idea to do this here. Isolate into its own method
    with open(exp_dir / 'compound_alternative_mass_dict.json', 'r') as jfile:
        side_product_dict = json.load(jfile)
    IS_compound_number = int(side_product_dict['IS_formula'].strip('SumF'))
    df_dict['IS'] = df_dict[IS_compound_number]
    del df_dict[IS_compound_number]

    return df_dict


def calculate_yield(dict_df):
    """calculate the mass response ratio (i.e. our approximation of yield) for a compounds vs. internal standard """
    # divide all values by the IS value
    for key, df in dict_df.items():
        df["yield"] = df["yield"].astype('float64') / dict_df["IS"]["yield"]
        df.fillna(value=0.0, inplace=True)  # TODO this might be a bad idea bc it hides errors. Better propagate na?

    return dict_df


def save_to_db(dict_df, db_path):
    # TODO this has too many functions. Isolate addition part.
    con = sqlite3.connect(db_path)
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


def import_lcms_unprocessed_data():
    df = pd.DataFrame
    return df


def calculate_lcms_yields(db_path, exp_dir):
    """
    Main function. Import lcms raw results from DB, apply evaluation logic (e.g. assign lcms peaks to products), and
    write results to DB table lcms_yields
    """
    #
    results_df = import_lcms_unprocessed_data()
    results_dict = split_products_into_dataframes(results_df, exp_dir)

    # calculate yields from areas
    yields = calculate_yield(results_dict)

    # save yields
    save_to_db(yields, db_path)


if __name__ == '__main__':
    calculate_lcms_yields(DB_PATH, EXP_DIR)
