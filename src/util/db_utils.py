"""
This module provides a class MyDatabaseConnection to facilitate querying the SQLite DB used for this project.
"""

import sqlite3 as sql
from typing import Optional, List, Tuple, Any

from PIL import Image
from rdkit.Chem import MolFromSmiles, Mol
from rdkit.Chem.Descriptors import MolWt

from src.definitions import DB_PATH


class MyDatabaseConnection:
    """
    This class provides methods to query the SQLite database.
    Popular queries are provided as methods, but arbitrary SELECT queries can be executed using the
    execute_arbitrary_simple_query.
    """

    def __init__(self):
        self.con = sql.connect(DB_PATH)
        self.cur = self.con.cursor()

    def execute_arbitrary_simple_query(self, query: str) -> list:
        """Execute simple SELECT statement on the database."""
        if "SELECT" not in query:
            raise ValueError("Not a SELECT statement")
        disallowed = [
            "INSERT",
            "CREATE",
            "ATTACH",
            "DETACH",
            "DROP",
            "UPDATE",
            "DELETE",
        ]
        for disallowed_statement in disallowed:
            if disallowed_statement.lower() in query.lower():
                raise ValueError("forbidden")
        return self.cur.execute(query).fetchall()

    def get_long_name(self, short: str) -> str:
        """Get long name (e.g. 2-Pyr001) from short name (e.g. I1))"""
        return self.cur.execute(
            "SELECT long FROM main.buildingblocks WHERE short = ?;", (short,)
        ).fetchone()[0]

    def get_smiles(self, short: str) -> str:
        """Get SMILES from a building block short"""
        smiles = self.cur.execute(
            "SELECT SMILES FROM main.buildingblocks WHERE short = ?;", (short,)
        ).fetchone()[0]
        return smiles

    def get_mol(self, short: str) -> Mol:
        """Get MOL from a building block short"""
        return MolFromSmiles(self.get_smiles(short))

    def show_image(self, short: str) -> None:
        """Show the molecular structure drawing for a building block short"""
        img_path = self.cur.execute(
            "SELECT image FROM main.buildingblocks WHERE short = ?;", (short,)
        ).fetchone()[0]
        Image.open(img_path).show()
        return

    def list_pg(self, short: str) -> tuple:
        """Returns a 4-tuple: (#boc, #cbz, #tbu, #tms)"""
        return self.cur.execute(
            "SELECT boc, cbz, tbu, tms FROM main.buildingblocks WHERE short = ?;",
            (short,),
        ).fetchone()

    def get_reactant_class(self, short: str) -> str:
        """Takes a building block short name and returns the reactant class"""
        return self.cur.execute(
            "SELECT reactant_class FROM main.buildingblocks WHERE short = ?;", (short,)
        ).fetchone()[0]

    def get_molecular_weight(self, short: str) -> float:
        """Get molecular weight from a building block short"""
        return MolWt(self.get_mol(short))

    def get_vl_member(self, vl_id: int) -> Mol:
        """Takes a vl_id and returns the product MOL"""
        smiles = self.cur.execute(
            "SELECT SMILES FROM main.virtuallibrary WHERE id = ?;", (vl_id,)
        ).fetchone()[0]
        return MolFromSmiles(smiles)

    def get_starting_materials_for_experiment(self, **kwargs: dict) -> tuple:
        """
        Return unique building blocks in an experiment. The experiment can be defined by either an exp_nr
        or a lab_journal_number.

        Args:
            kwargs (dict): Exactly one key of lab_journal_nr or exp_nr must be given.

        Returns:
            tuple, 3-tuple of sorted building block lists: ([initiators], [monomers], [terminators])
        """
        # if-clause checks whether exactly one of the possible two kwargs(and no other kwargs) was given
        if not set() < kwargs.keys() < {"lab_journal_number", "exp_nr"}:
            raise ValueError(
                "One keyword argument is required: lab_journal_number=x or exp_nr=x"
            )

        if list(kwargs.keys())[0] == "exp_nr":
            initiators = sorted(
                [
                    i[0]
                    for i in self.cur.execute(
                        f"SELECT DISTINCT initiator FROM experiments WHERE exp_nr = ?;",
                        (kwargs["exp_nr"],),
                    ).fetchall()
                ],
                key=lambda x: int(x[1:]),
            )
            monomers = sorted(
                [
                    i[0]
                    for i in self.cur.execute(
                        f"SELECT DISTINCT monomer FROM experiments WHERE exp_nr = ?;",
                        (kwargs["exp_nr"],),
                    ).fetchall()
                ],
                key=lambda x: int(x[1:]),
            )
            terminators = sorted(
                [
                    i[0]
                    for i in self.cur.execute(
                        f"SELECT DISTINCT terminator FROM experiments WHERE exp_nr = ?;",
                        (kwargs["exp_nr"],),
                    ).fetchall()
                ],
                key=lambda x: int(x[1:]),
            )
        else:
            initiators = sorted(
                [
                    i[0]
                    for i in self.cur.execute(
                        f"SELECT DISTINCT initiator FROM experiments WHERE lab_journal_number = ?;",
                        (kwargs["lab_journal_number"],),
                    ).fetchall()
                ],
                key=lambda x: int(x[1:]),
            )
            monomers = sorted(
                [
                    i[0]
                    for i in self.cur.execute(
                        f"SELECT DISTINCT monomer FROM experiments WHERE lab_journal_number = ?;",
                        (kwargs["lab_journal_number"],),
                    ).fetchall()
                ],
                key=lambda x: int(x[1:]),
            )
            terminators = sorted(
                [
                    i[0]
                    for i in self.cur.execute(
                        f"SELECT DISTINCT terminator FROM experiments WHERE lab_journal_number = ?;",
                        (kwargs["lab_journal_number"],),
                    ).fetchall()
                ],
                key=lambda x: int(x[1:]),
            )
        return initiators, monomers, terminators

    def get_product_mol_by_well(
        self, lab_journal_number: str, well: str, product_type: str
    ) -> Optional[Mol]:
        """
        Args:
            lab_journal_number (str): Unique identifier for the plate
            well (str): Well identifier within the plate
            product_type (str): [A-H]

        Returns:
            Mol (optional): rdkit mol object, or None if no match was found
        """
        product_mapping = dict(zip("ABCDEFGH", range(8)))
        product_idx = product_mapping[product_type]
        result = self.cur.execute(
            "SELECT product_A_smiles, product_B_smiles, product_C_smiles, product_D_smiles, product_E_smiles, \
            product_F_smiles, product_G_smiles, product_H_smiles \
            FROM experiments WHERE lab_journal_number = ? AND well = ?;",
            (lab_journal_number, well),
        ).fetchall()
        try:
            smiles = result[0][product_idx]
        except KeyError:
            return None
        return MolFromSmiles(smiles)

    def get_number_of_experiments_by_date(self) -> List[Tuple[str, int]]:
        return self.cur.execute(
            "SELECT synthesis_date_unixepoch, COUNT(*) FROM experiments GROUP BY synthesis_date_unixepoch;"
        ).fetchall()

    def get_experiments_with_buildingblock(self, short) -> List[Tuple[Any]]:
        """
        Returns a list of experiments that contain the given building block.
        Each entry is a tuple consisting of exp_nr, plate_nr, well, lab_journal_number, initiator, monomer, terminator,
        product_A_lcms_ratio, product_B_lcms_ratio, product_C_lcms_ratio, product_D_lcms_ratio, product_E_lcms_ratio,
        product_F_lcms_ratio, product_G_lcms_ratio, product_H_lcms_ratio, vl_id, valid
        """
        return self.cur.execute(
            "SELECT exp_nr, plate_nr, well, lab_journal_number, initiator, monomer, terminator, product_A_lcms_ratio, product_B_lcms_ratio, product_C_lcms_ratio, product_D_lcms_ratio, product_E_lcms_ratio, product_F_lcms_ratio, product_G_lcms_ratio, product_H_lcms_ratio, vl_id, valid FROM experiments WHERE initiator = ? OR monomer = ? OR terminator = ?;",
            (short, short, short),
        ).fetchall()

    def __delete__(self, instance):
        self.con.close()
