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
import pandas as pd
from config import *

exp_nr = 'JG221'
exp_dir = PLATES_DIR / exp_nr


def import_lcms_unprocessed_data(db_cur, exp_nr):
    # TODO this is a terrible but working solution
    # for now, we do this for a certain synthesis id. we could also do it for the entire DB
    res = db_cur.execute(
        'SELECT synthesis_id, lcms_compounds, lcms_areas FROM lcms WHERE synthesis_id IN (SELECT id FROM experiments WHERE lab_journal_number = ?);',
        (exp_nr,)).fetchall()
    # for now, we produce a df to exchange our data. Other formats might be more efficient.
    first_res = [
        [int(res[0][0])] + [float(i.strip('[').strip(']')) for i in res[0][2].split(',')]]  # TODO horrible parsing
    first_header = ['synthesis_id'] + [i.split("'")[1] for i in res[0][1].split(
        ',')]  # TODO this is a ridiculousy breakable way to parse that string
    df = pd.DataFrame(data=first_res, columns=first_header)
    # now that we have the df, we iterate all further results.
    for result in res[1:]:
        columns = ['synthesis_id'] + [i.split("'")[1] for i in result[1].split(',')]  # again, the parsing is horrible
        values = [[int(result[0])] + [float(i.strip('[').strip(']')) for i in result[2].split(',')]]
        df_temp = pd.DataFrame(data=values, columns=columns)
        df = df.append(df_temp)
        # for c,v in zip(columns, values)
    # we need to get the plate, row, column values for the experiment db
    synthesis_ids = df['synthesis_id'].tolist()
    res = db_cur.execute(
        f'SELECT id, plate_nr, well FROM main.experiments WHERE id IN ({",".join(["?"] * len(synthesis_ids))})',
        synthesis_ids).fetchall()

    res = [[r[0], r[1], r[2][0], r[2][1:]] for r in res]
    # rows = [w[0] for w in wells]
    # columns = [w[1:] for w in wells]
    # entries = [[i, p, r, c] for i,p,r,c in zip(ids, plates, rows, columns)]
    df_well = pd.DataFrame(columns=['synthesis_id', 'plate', 'row', 'column'], data=res)
    df = df.join(df_well.set_index('synthesis_id'), on='synthesis_id', how='inner')
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
        df_dict[number] = df[['plate', 'row', 'column', c]]
        df_dict[number].rename(columns={c: 'yield'}, inplace=True)

    # Add df for BPC to dict as well
    df_dict['BPC'] = df[['plate', 'row', 'column', 'BPC Area']].rename(columns={'BPC Area': 'yield'})

    # Generate dataframe for IS TODO bad idea to do this here. Isolate into its own method
    with open(exp_dir / 'compound_alternative_mass_dict.json', 'r') as jfile:
        side_product_dict = json.load(jfile)
    # TODO remove this as soon as we are only using PG stuff
    # side_product_dict = {"A_formula": "SumF1", "B_formula": "SumF2", "C_formula": "SumF3", "D_formula": "SumF4", "E_formula": "SumF5", "F_formula": "SumF6", "G_formula": "SumF7", "H_formula": "SumF8", "IS_formula": "SumF9"}

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


def save_to_db(dict_df, db_path, exp_nr):
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
        for i, data in product_yields.iterrows():  # TODO should not write to exp db but have its own lcms yields thing
            cur.execute(
                f'UPDATE main.experiments SET product_{k}_lcms_ratio = ? WHERE lab_journal_number = ? AND well = ?;',
                (data['yield'], exp_nr, f'{data["row"]}{data["column"]}'))
        con.commit()
    con.close()
    return


def calculate_lcms_yields(db_path, exp_dir, exp_nr):
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
    results_df = import_lcms_unprocessed_data(cur, exp_nr)
    results_dict = split_products_into_dataframes(results_df, exp_dir)

    # calculate yields from areas
    yields = calculate_yield(results_dict)

    # save yields
    save_to_db(yields, db_path, exp_nr)
    con.close()


if __name__ == '__main__':
    calculate_lcms_yields(DB_PATH, exp_dir, exp_nr)
