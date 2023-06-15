"""
This module provides a class SynFermDatabaseConnection to facilitate querying the SQLite DB used for this project.
"""

import sqlite3 as sql
from typing import Optional, List, Tuple, Any, Union

from PIL import Image
from rdkit.Chem import MolFromSmiles, Mol
from rdkit.Chem.Descriptors import MolWt
import pandas as pd

from src.definitions import DB_PATH


class SynFermDatabaseConnection:
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

    def get_reaction_id(
        self, identifier: Union[int, Tuple[str, str], Tuple[int, int, str]]
    ) -> int:
        """
        Get reaction id (primary key in experiments table) from a tuple of (lab_journal_number, well)
        or (exp_nr, plate_nr, well). If an integer is passed, checks if it is a valid reaction id.

        Args:
            identifier: tuple of (lab_journal_number, well) or (exp_nr, plate_nr, well) or a single integer (reaction id)

        Returns:
            int: ID of the reaction instance in the experiments table

        Raises:
            ValueError: if an integer is passed that is not a valid reaction id
        """
        if isinstance(identifier, int):
            if self.cur.execute(
                "SELECT COUNT(*) FROM main.experiments WHERE id = ?;", (identifier,)
            ).fetchone()[0]:
                return identifier
            else:
                raise ValueError("Reaction ID not found in database")

        if len(identifier) == 2:
            return self.cur.execute(
                "SELECT id FROM main.experiments WHERE lab_journal_number = ? AND well = ?;",
                identifier,
            ).fetchone()[0]
        elif len(identifier) == 3:
            return self.cur.execute(
                "SELECT id FROM main.experiments WHERE exp_nr = ? AND plate_nr = ? AND well = ?;",
                identifier,
            ).fetchone()[0]
        else:
            raise ValueError("Identifier must be a tuple of length 2 or 3")

    def get_long_name(self, short: str) -> str:
        """Get long name (e.g. 2-Pyr001) from short name (e.g. I1))"""
        return self.cur.execute(
            "SELECT long FROM main.buildingblocks WHERE short = ?;", (short,)
        ).fetchone()[0]

    def get_reaction_ids_for_plate(
        self, identifier: Union[str, Tuple[int, int]]
    ) -> List[int]:
        """
        Get reaction ids (primary keys in experiments table) from a string (lab_journal_number)
        or a 2-tuple (exp_nr, plate_nr).

        Args:
            identifier: string (lab_journal_number) or tuple (exp_nr, plate_nr)

        Returns:
            List[int]: IDs of the reaction instances connected to the given identifier in the experiments table

        Raises:
            ValueError: If an invalid identifier is passed
        """

        if isinstance(identifier, str):
            return [
                x[0]
                for x in self.cur.execute(
                    "SELECT id FROM main.experiments WHERE lab_journal_number = ?;",
                    (identifier,),
                ).fetchall()
            ]
        elif isinstance(identifier, tuple):
            return [
                x[0]
                for x in self.cur.execute(
                    "SELECT id FROM main.experiments WHERE exp_nr = ? AND plate_nr = ?;",
                    identifier,
                ).fetchall()
            ]
        else:
            raise ValueError("Identifier must be a string or a tuple")

    def get_reaction_ids_for_building_block(
        self,
        initiator: str = "I%",
        monomer: str = "M%",
        terminator: str = "T%",
        filter_exp_nr: Tuple[int, int] = (1, 99999),
    ) -> List[int]:
        """
        Get reaction ids (primary keys in experiments table) for a building block, or combination of building blocks.
        If no building block is specified, all reaction ids are returned.

        Args:
            initiator (str): short name of the initiator building block
            monomer (str): short name of the monomer building block
            terminator (str): short name of the terminator building block

        Returns:
            List[int]: IDs of the reaction instances connected to the given building blocks in the experiments table

        Raises:
            ValueError: If no building block is specified
        """
        return [
            x[0]
            for x in self.cur.execute(
                "SELECT id FROM experiments WHERE initiator LIKE ? AND monomer LIKE ? AND terminator LIKE ? AND exp_nr BETWEEN ? AND ?;",
                (initiator, monomer, terminator, filter_exp_nr[0], filter_exp_nr[1]),
            ).fetchall()
        ]

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

    def get_building_block_class(self, short: str) -> str:
        """
        Takes a building block short name and returns the building block class.

        Args:
            short (str): building block short name

        Returns:
            str: building block class. If short is not found, returns None
        """
        res = self.cur.execute(
            "SELECT reactant_class FROM main.buildingblocks WHERE short = ?;", (short,)
        ).fetchone()
        try:
            return res[0]
        except TypeError:
            return None

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

    def get_starting_materials_for_reaction(self, identifier) -> list:
        """
        Return the SMILES of the starting materials for a given reaction.

        Args:
            identifier (int): Reaction identifier. This can be either
                a single int (id in the experiments table) or
                a 2-tuple of str (lab_journal_nr, well) or
                a 3-tuple of two int and one str (exp_nr, plate_nr, well)

        Returns:
            list: list of starting material SMILES ([initiator, monomer, terminator])
        """
        reaction_id = self.get_reaction_id(identifier)
        result = self.cur.execute(
            """
            SELECT initiator, monomer, terminator
            FROM experiments
            WHERE id = ?;
            """,
            (reaction_id,),
        ).fetchone()

        if result is None:
            raise ValueError(f"No reaction found for {identifier}")

        return [self.get_smiles(i) for i in result]

    def get_product_mol(
        self,
        identifier: Union[int, Tuple[str, str], Tuple[int, int, str]],
        product_type: str,
    ) -> Optional[Mol]:
        """
        Get the product mol for a given well in a given plate.
        Args:
            identifier (int or tuple): Either the reaction id or a tuple of (lab_journal_number, well)
                or (exp_nr, plate_nr, well)
            product_type (str): [A-H]

        Returns:
            Mol (optional): rdkit mol object, or None if no match was found
        """
        product_mapping = dict(zip("ABCDEFGH", range(8)))
        result = self.get_product_smiles(identifier)

        for pt in product_type:
            try:
                smiles = result[product_mapping[pt]]
            except KeyError:
                yield None
            yield MolFromSmiles(smiles)

    def get_product_smiles(
        self, identifier: Union[int, Tuple[str, str], Tuple[int, int, str]]
    ) -> Tuple[str]:
        """
        Get all product formulae for a given experiment and well.
        Args:
            identifier (int or tuple): Either the reaction id or a tuple of (lab_journal_number, well)
                or (exp_nr, plate_nr, well)

        Returns:
            tuple of product formulae
        """
        reaction_id = self.get_reaction_id(identifier)

        return self.cur.execute(
            "SELECT product_A_smiles, product_B_smiles, product_C_smiles, product_D_smiles, product_E_smiles, \
            product_F_smiles, product_G_smiles, product_H_smiles \
            FROM experiments WHERE id = ?;",
            (reaction_id,),
        ).fetchall()[0]

    def get_lcms_peaks(
        self,
        identifier,
        with_delta=False,
        with_assignment=False,
        with_building_blocks=False,
    ) -> pd.DataFrame:
        """
        Get the LCMS peaks (as extracted from PDF) for a given reaction.

        Args:
            identifier (int or tuple): Either the reaction id or a tuple of (lab_journal_number, well) or (exp_nr, plate_nr, well)
            with_delta (bool): If True, also return the mass differences for different I, M, T combinations for each peak.
                In this case, only peaks for which the difference was calculated are returned.
            with_assignment (bool): If True, also return the assignments for each peak, if available.

        Returns:
            pd.DataFrame: DataFrame of the peaks, optionally with the mass differences
        """
        reaction_id = self.get_reaction_id(identifier)
        columns = [
            "lcms_peaks.id",
            "lcms_peaks.experiment_id",
            "lcms_peaks.peak_nr",
            "retention_time_s",
            "area",
            "intensity",
            "signal_to_noise",
            "mz_max",
            "fwhm_min",
            '"%area"',
            '"%intensity"',
        ]
        if with_building_blocks:
            columns.extend(["initiator", "monomer", "terminator"])
        if with_delta:
            columns.extend(
                [
                    "d.delta_I",
                    "d.delta_M",
                    "d.delta_T",
                    "d.delta_Iacid",
                    "d.delta_bAA",
                    "d.delta_IM",
                    "d.delta_IT",
                    "d.delta_MT",
                    "d.delta_IMT",
                ]
            )
        if with_assignment:
            columns.append("a.assignment")

        query = f"SELECT {', '.join(columns)} FROM lcms_peaks"
        if with_building_blocks:
            query += " LEFT JOIN experiments e on lcms_peaks.experiment_id = e.id"
        if with_delta:
            query += " JOIN lcms_peaks_differences d on lcms_peaks.id = d.peak_id"
        if with_assignment:
            query += " LEFT JOIN lcms_peaks_assignment a on lcms_peaks.id = a.peak_id"
        query += " WHERE lcms_peaks.experiment_id = ?;"
        result = self.cur.execute(query, (reaction_id,)).fetchall()
        df = pd.DataFrame(result, columns=[col.split(".")[-1] for col in columns])
        return df

    def get_experiments_table_as_df(self) -> pd.DataFrame:
        """Returns the experiments table as a pandas.Dataframe"""
        return pd.read_sql_query(
            "SELECT * FROM experiments", self.con, index_col="id", coerce_float=False
        )

    def get_number_of_experiments_by_date(self) -> List[Tuple[str, int]]:
        """Returns a list of tuples consisting of (date, number of experiments on that date)"""
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
