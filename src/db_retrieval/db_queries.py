"""
Module to query our specific DB implementation.
Provides only functions for import, to facilitate querying the DB
"""

import sqlite3 as sql

from rdkit import Chem
from PIL import Image

from definitions import DB_PATH


class MyDatabaseConnection:

    def __init__(self):
        self.con = sql.connect(DB_PATH)
        self.cur = self.con.cursor()

    """Some functions for simple queries to the db"""

    def get_mol(self, short):
        smiles = self.cur.execute('SELECT SMILES FROM main.buildingblocks WHERE short = ?;', (short,)).fetchone()[0]
        return Chem.MolFromSmiles(smiles)

    def get_smiles(self, short):
        smiles = self.cur.execute('SELECT SMILES FROM main.buildingblocks WHERE short = ?;', (short,)).fetchone()[0]
        return smiles

    def show_image(self, short):
        img_path = self.cur.execute('SELECT image FROM main.buildingblocks WHERE short = ?;', (short,)).fetchone()[0]
        Image.open(img_path).show()
        return

    def list_pg(self, short):
        """Returns a 4-tuple: (#boc, #cbz, #tbu, #tms)"""
        return self.cur.execute('SELECT boc, cbz, tbu, tms FROM main.buildingblocks WHERE short = ?;',
                                (short,)).fetchone()

    def __delete__(self, instance):
        self.con.close()


if __name__ == '__main__':
    # add debugging statements here
    pass
