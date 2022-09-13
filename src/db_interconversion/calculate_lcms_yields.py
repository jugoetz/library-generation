"""
Go from raw lcms analysis results (i.e. peak areas) to a yield metric (i.e. peak area / IS area for every compound).

We start with a DB table 'lcms':
id | synthesis id  | [lcms compounds] | [lcms areas]

We add data into the DB table 'experiments':
id (foreign key from lcms table) | SumF.. | type (e.g. A) | SMILES (or better: id of the product in a different table)

Note that this overwrites previous data on every execution

This script raises a RuntimeError if there are multiple entries for a given synthesis_id in the lcms table.

edit config.yaml before running
"""

import re
import sqlite3
from functools import reduce

import pandas as pd

from definitions import PLATES_DIR, DB_PATH
from utils import get_product_dict, get_internal_standard_number, get_conf

# configurations
# edit config.yaml before running
conf = get_conf()


def strip_brackets(s):
    """Strip brackets from a string."""
    return s.strip('[').strip(']')


def import_lcms_for_plate(db_cur, lab_journal_nr, unique_ids=True):
    """
    Import the lcms data for a given plate from the DB. Return a dataframe with the data.
    Since this is meant for one plate, we require that the positional meaning of the data is consistent.
    Returns a df indexed by synthesis_id, containing:
     - an arbitrary number of data columns
     - one column holding the plate number
     - one column holding the row letter
     - one column holding the column number
     If unique_ids is True, a RuntimeError is raised if there are multiple entries for a given synthesis_id.
    """
    # for now, we do this for a certain synthesis id. we could also do it for the entire DB
    res = db_cur.execute(
        'SELECT synthesis_id, lcms_compounds, lcms_areas FROM lcms WHERE synthesis_id IN (SELECT id FROM experiments WHERE lab_journal_number = ?);',
        (lab_journal_nr,)).fetchall()
    # we rearrange the data to put it into a dataframe
    synthesis_ids, headers, data = [[x[i] for x in res] for i in range(3)]
    # check that all headers are identical
    if not headers.count(headers[0]) == len(headers):
        raise ValueError('The data headers are not identical for all synthesis ids')
    # check that all synthesis ids are unique
    if unique_ids and not len(synthesis_ids) == len(set(synthesis_ids)):
        raise RuntimeError('There are multiple entries for a given synthesis id')
    # cast data to float
    data = [[float(strip_brackets(s)) for s in d.split(",")] for d in data]
    # put everything into a dataframe
    df = pd.DataFrame(data=data, columns=eval(headers[0]), index=synthesis_ids)
    return df


def split_products_into_dataframes(df, exp_dir):
    """Generate individual dataframe for each product (and for special cases like IS and BPC). Return a dict of those"""
    df_dict = {}

    # find the product columns in df by regex
    regex = re.compile('SumF([0-9]+) Area')
    product_columns = filter(regex.match, df.columns)

    # for every product column, add a separate df to dict
    for c in product_columns:
        number = int(regex.match(c).group(1))
        df_dict[number] = df.loc[:, [c]]
        df_dict[number] = df_dict[number].rename(columns={c: 'yield'})

    # Add df for BPC to dict as well
    df_dict['BPC'] = df.loc[:, ['BPC Area']].rename(columns={'BPC Area': 'yield'})

    # Generate dataframe for IS
    internal_standard_number = int(get_internal_standard_number(exp_dir).strip('SumF'))
    df_dict['IS'] = df_dict[internal_standard_number]
    del df_dict[internal_standard_number]

    return df_dict


def calculate_lcms_yield(dict_df):
    """Calculate the mass response ratio (i.e. our approximation of yield) for a compounds vs. internal standard """
    # divide all values by the IS value
    for key, df in dict_df.items():
        df.loc[:, "yield"] = df.loc[:, "yield"].astype('float64') / dict_df["IS"].loc[:, "yield"]

    return dict_df


def save_lcms_yields_to_db(dict_df, db_path):
    """From the LCMS yield data, extract the yield for each product and save it to the DB"""
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    # load the saved connections between product type (e.g. A) and SumFxx field
    product_dict = get_product_dict(exp_dir)
    # reformat dict for easier use
    product_associations = {t: [v for k, v in product_dict.items() if k.startswith(t)] for t in 'ABCDEFGH'}
    # sum over all dataframes for a product type
    for k, v in product_associations.items():
        # get all dataframes for a product type
        product_yields = [dict_df[int(val.strip('SumF'))] for val in v]
        # sum over all dataframes for a product type
        yield_df = reduce(lambda a, b: a.add(b, fill_value=0.0), product_yields)
        # save to DB
        for i, data in yield_df.iterrows():
            cur.execute(
                f'UPDATE main.experiments SET product_{k}_lcms_ratio = ? WHERE id = ?;',
                (data['yield'], i))
        con.commit()
    con.close()
    return


def calculate_lcms_yields(db_path, exp_dir, lab_journal_nr):
    """
    Main function. Import lcms raw results from DB, apply evaluation logic (e.g. assign lcms peaks to products), and
    write results to DB table lcms_yields
    What's part of the evaluation logic:
    - Identify the IS peak area
    - Identify and discard BPC area
    - Normalize all other peaks by diving by IS peak area
    - Map LCMS output names (e.g. SumF1 to the product that is detected. This should be identified with an experiment product (e.g. Prod F - 2 Boc) AND an actual SMILES)
    """
    con = sqlite3.connect(DB_PATH)
    cur = con.cursor()
    results_df = import_lcms_for_plate(cur, lab_journal_nr)
    results_dict = split_products_into_dataframes(results_df, exp_dir)

    # calculate yields from areas
    yields = calculate_lcms_yield(results_dict)

    # save yields
    save_lcms_yields_to_db(yields, db_path)
    con.close()


if __name__ == '__main__':
    for lab_journal_nr in conf['lab_journal_numbers']:
        print(f'Now calculating LCMS ratios for {lab_journal_nr}...')
        exp_dir = PLATES_DIR / lab_journal_nr
        calculate_lcms_yields(DB_PATH, exp_dir, lab_journal_nr)
    print('Calculation of LCMS yields finished!')
