"""
Module to query our specific DB implementation.
Provides only functions for import, to facilitate querying the DB
"""

import sqlite3 as sql

from rdkit.Chem import MolFromSmiles
from PIL import Image

from definitions import DB_PATH


class MyDatabaseConnection:

    def __init__(self):
        self.con = sql.connect(DB_PATH)
        self.cur = self.con.cursor()

    """Some functions for simple queries to the db"""

    def get_smiles(self, short: str):
        """Get SMILES from a building block short"""
        smiles = self.cur.execute('SELECT SMILES FROM main.buildingblocks WHERE short = ?;', (short,)).fetchone()[0]
        return smiles

    def get_mol(self, short: str):
        """Get MOL from a building block short"""
        return MolFromSmiles(self.get_smiles(short))

    def show_image(self, short: str):
        """Show the molecular structure drawing for a building block short"""
        img_path = self.cur.execute('SELECT image FROM main.buildingblocks WHERE short = ?;', (short,)).fetchone()[0]
        Image.open(img_path).show()
        return

    def list_pg(self, short: str):
        """Returns a 4-tuple: (#boc, #cbz, #tbu, #tms)"""
        return self.cur.execute('SELECT boc, cbz, tbu, tms FROM main.buildingblocks WHERE short = ?;',
                                (short,)).fetchone()

    def get_reactant_class(self, short: str):
        """Takes a building block short name and returns the reactant class"""
        return self.cur.execute('SELECT reactant_class FROM main.buildingblocks WHERE short = ?;', (short,)).fetchone()[
            0]

    def get_vl_member(self, vl_id):
        """Takes a vl_id and returns the product MOL"""
        smiles = self.cur.execute('SELECT SMILES FROM main.virtuallibrary WHERE id = ?;',
                                  (vl_id,)).fetchone()[0]
        return MolFromSmiles(smiles)

    def get_product_mol_by_well(self, lab_journal_number: str, well: str, product_type: str):
        """
        :param lab_journal_number: Unique identifier for the plate
        :param well: Identifier for the well within the plate
        :param product_type: [A-H]
        :return: rdkit mol object, or None if no match was found
        """
        product_mapping = dict(zip('ABCDEFGH', range(8)))
        product_idx = product_mapping[product_type]
        result = self.cur.execute(
            'SELECT product_A_smiles, product_B_smiles, product_C_smiles, product_D_smiles, product_E_smiles, \
            product_F_smiles, product_G_smiles, product_H_smiles \
            FROM experiments WHERE lab_journal_number = ? AND well = ?;',
            (lab_journal_number, well)).fetchall()
        try:
            smiles = result[0][product_idx]
        except KeyError:
            return None
        return MolFromSmiles(smiles)

    def __delete__(self, instance):
        self.con.close()


if __name__ == '__main__':
    # add debugging statements here
    pass
