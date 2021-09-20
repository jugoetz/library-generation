import sqlite3

from rdkit import Chem
from rdkit.Chem import Draw

from definitions import DB_PATH


def get_product(lab_journal_number, well, product_type):
    product_mapping = dict(zip('ABCDEFGH', range(1, 9)))
    product_idx = product_mapping[product_type]
    # query db
    con = sqlite3.connect(DB_PATH)
    cur = con.cursor()
    result = cur.execute(
        'SELECT long_name, product_A_smiles, product_B_smiles, product_C_smiles, product_D_smiles, product_E_smiles, product_F_smiles, product_G_smiles, product_H_smiles FROM experiments WHERE lab_journal_number = ? AND well = ?;',
        (lab_journal_number, well)).fetchall()
    if len(result) == 0:
        print('Error: No entry in DB for the given parameters')
    long_name = result[0][0]
    smiles = result[0][product_idx]
    print(long_name)
    print(smiles)
    if smiles is not None:
        Draw.MolToImage(Chem.MolFromSmiles(smiles), size=(500, 500)).show()
    else:
        print('This molecule does not exist')

    con.close()
    return


if __name__ == '__main__':
    lab_journal_number = 'JG264'
    well = 'A5'
    product = 'A'
    get_product(lab_journal_number, well, product)
