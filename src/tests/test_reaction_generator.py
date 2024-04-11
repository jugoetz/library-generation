import unittest
import random

from rdkit import Chem, RDLogger

RDLogger.DisableLog(
    "rdApp.warning"
)  # suppress warnings about unmapped atoms in reactionSMILES

from src.library_design.reaction_generator import SFReactionGenerator
from src.util.db_utils import SynFermDatabaseConnection
from src.util.rdkit_util import desalt_building_block, remove_monomer_pg_chirality
from src.definitions import DB_PATH

DB_EXISTS = DB_PATH.exists()


class TestSFReactionGenerator(unittest.TestCase):
    """
    Note that most of these tests rely on the database being present and populated. They will be skipped if the database is not present.
    """

    def setUp(self) -> None:
        self.rxn_generator = SFReactionGenerator()
        self.con = SynFermDatabaseConnection()
        self.spiro_c_pattern = Chem.MolFromSmarts(
            "[$([CR2](O1)(ONC2)(C2)C(=O)OC1)]"
        )  # hits a spiro carbon atom between an isoxazolidine and a 5-membered lactone-ether ring
        if DB_EXISTS:
            self.building_blocks = {
                long: Chem.MolToSmiles(
                    remove_monomer_pg_chirality(desalt_building_block(smiles))
                )
                if (
                    long.startswith("Mon")
                    or long.startswith("Fused")
                    or long.startswith("Spiro")
                )
                else Chem.MolToSmiles(desalt_building_block(smiles))
                for long, smiles in self.con.building_blocks()
            }

    @unittest.skipIf(not DB_EXISTS, "Missing database")
    def test_gives_correct_starting_materials(self):
        """
        Test that the reaction generator gives the correct starting materials for a given product
        on a large random subset of the VL.
        """
        res = self.con.con.execute(
            "SELECT id, initiator_long, monomer_long, terminator_long, SMILES FROM virtuallibrary WHERE type = 'A' AND initiator_long != '4-Pyrazole002';",
        ).fetchall()  # 4-Pyrazole002 is the bullshit initiator with "urea-KAT"

        for row in random.sample(res, 10000):  # we test 10000 random VL members
            try:
                generated_reactants = self.rxn_generator.generate_reactants(row[4])
            except (ValueError, RuntimeError) as e:
                with self.subTest(vl_id=row[0]):
                    self.fail(f"Error for vl_id {row[0]}: {e}")

            with self.subTest(bb=row[1], vl_id=row[0]):
                self.assertEqual(
                    self.building_blocks[row[1]],
                    Chem.MolToSmiles(generated_reactants[0]),
                )
            with self.subTest(bb=row[2], vl_id=row[0]):
                self.assertEqual(
                    self.building_blocks[row[2]],
                    Chem.MolToSmiles(generated_reactants[1]),
                )
            with self.subTest(bb=row[3], vl_id=row[0]):
                self.assertEqual(
                    self.building_blocks[row[3]],
                    Chem.MolToSmiles(generated_reactants[2]),
                )

    @unittest.skipIf(not DB_EXISTS, "Missing database")
    def test_gives_correct_product_A(self):
        """
        Test that the reaction generator gives the correct product for a given set of starting materials
        on a large random subset of the VL.
        """
        res = self.con.con.execute(
            "SELECT id, initiator_long, monomer_long, terminator_long, SMILES FROM virtuallibrary WHERE type = 'A' AND initiator_long != '4-Pyrazole002';",
        ).fetchall()  # 4-Pyrazole002 is the bullshit initiator with "urea-KAT"

        for row in random.sample(res, 10000):  # we test 10000 random VL members
            try:
                generated_product = Chem.MolToSmiles(
                    self.rxn_generator.generate_product(
                        [
                            Chem.MolFromSmiles(self.building_blocks[bb])
                            for bb in row[1:4]
                        ]
                    )
                )
            except (ValueError, RuntimeError) as e:
                with self.subTest(vl_id=row[0], msg="Exception raised"):
                    self.fail(f"Error for vl_id {row[0]}: {e}")

            with self.subTest(vl_id=row[0]):
                self.assertEqual(
                    row[4],
                    generated_product,
                )

    @unittest.skipIf(not DB_EXISTS, "Missing database")
    def test_gives_correct_product_B(self):
        """
        Test that the reaction generator gives the correct product for a given set of starting materials
        on a large random subset of the VL.
        """
        res = self.con.con.execute(
            "SELECT id, initiator_long, monomer_long, terminator_long, SMILES FROM virtuallibrary WHERE type = 'B' AND initiator_long != '4-Pyrazole002';",
        ).fetchall()  # 4-Pyrazole002 is the bullshit initiator with "urea-KAT"

        for row in random.sample(res, 10000):  # we test 10000 random VL members
            try:
                generated_product = Chem.MolToSmiles(
                    self.rxn_generator.generate_product(
                        [
                            Chem.MolFromSmiles(self.building_blocks[bb])
                            for bb in row[1:4]
                        ],
                        product_type="B",
                    )
                )
            except (ValueError, RuntimeError) as e:
                with self.subTest(vl_id=row[0], msg="Exception raised"):
                    self.fail(f"Error for vl_id {row[0]}: {e}")

            with self.subTest(vl_id=row[0]):
                self.assertEqual(
                    row[4],
                    generated_product,
                )

    @unittest.skipIf(not DB_EXISTS, "Missing database")
    def test_gives_correct_product_C(self):
        """
        Test that the reaction generator gives the correct product for a given set of starting materials
        on a large random subset of the VL.
        """
        res = self.con.con.execute(
            "SELECT id, initiator_long, monomer_long, terminator_long, SMILES FROM virtuallibrary WHERE type = 'C' AND initiator_long != '4-Pyrazole002';",
        ).fetchall()  # 4-Pyrazole002 is the bullshit initiator with "urea-KAT"

        for row in random.sample(res, 10000):  # we test 10000 random VL members
            try:
                generated_product = Chem.MolToSmiles(
                    self.rxn_generator.generate_product(
                        [
                            Chem.MolFromSmiles(self.building_blocks[bb])
                            for bb in row[1:4]
                        ],
                        product_type="C",
                    )
                )
            except (ValueError, RuntimeError) as e:
                with self.subTest(vl_id=row[0], msg="Exception raised"):
                    self.fail(f"Error for vl_id {row[0]}: {e}")

            with self.subTest(vl_id=row[0]):
                self.assertEqual(
                    row[4],
                    generated_product,
                )

    @unittest.skipIf(not DB_EXISTS, "Missing database")
    def test_gives_correct_product_D(self):
        """
        Test that the reaction generator gives the correct product for a given set of starting materials
        on a large random subset of the VL.
        """
        res = self.con.con.execute(
            "SELECT id, initiator_long, monomer_long, terminator_long, SMILES FROM virtuallibrary WHERE type = 'D' AND initiator_long != '4-Pyrazole002';",
        ).fetchall()  # 4-Pyrazole002 is the bullshit initiator with "urea-KAT"

        for row in random.sample(res, 10000):  # we test 10000 random VL members
            try:
                generated_product = Chem.MolToSmiles(
                    self.rxn_generator.generate_product(
                        [
                            Chem.MolFromSmiles(self.building_blocks[bb])
                            for bb in row[1:4]
                        ],
                        product_type="D",
                    )
                )
            except (ValueError, RuntimeError) as e:
                with self.subTest(vl_id=row[0], msg="Exception raised"):
                    self.fail(f"Error for vl_id {row[0]}: {e}")

            with self.subTest(vl_id=row[0]):
                self.assertEqual(
                    row[4],
                    generated_product,
                )

    @unittest.skipIf(not DB_EXISTS, "Missing database")
    def test_gives_correct_product_E(self):
        """
        Test that the reaction generator gives the correct product for a given set of starting materials
        on a large random subset of the VL.
        """
        res = self.con.con.execute(
            "SELECT id, initiator_long, monomer_long, terminator_long, SMILES FROM virtuallibrary WHERE type = 'E' AND initiator_long != '4-Pyrazole002';",
        ).fetchall()  # 4-Pyrazole002 is the bullshit initiator with "urea-KAT"

        for row in random.sample(res, 10000):  # we test 10000 random VL members
            try:
                generated_product = Chem.MolToSmiles(
                    self.rxn_generator.generate_product(
                        [
                            Chem.MolFromSmiles(self.building_blocks[bb])
                            for bb in row[1:4]
                        ],
                        product_type="E",
                    )
                )
            except (ValueError, RuntimeError) as e:
                with self.subTest(vl_id=row[0], msg="Exception raised"):
                    self.fail(f"Error for vl_id {row[0]}: {e}")

            with self.subTest(vl_id=row[0]):
                self.assertEqual(
                    row[4],
                    generated_product,
                )

    @unittest.skipIf(not DB_EXISTS, "Missing database")
    def test_gives_correct_product_F(self):
        """
        Test that the reaction generator gives the correct product for a given set of starting materials
        on a large random subset of the VL.
        """

        res = self.con.con.execute(
            "SELECT id, initiator_long, monomer_long, terminator_long, SMILES FROM virtuallibrary WHERE type = 'F' AND initiator_long != '4-Pyrazole002';",
        ).fetchall()  # 4-Pyrazole002 is the bullshit initiator with "urea-KAT"

        for row in random.sample(res, 10000):  # we test 10000 random VL members
            try:
                generated_product = Chem.MolToSmiles(
                    self.rxn_generator.generate_product(
                        [
                            Chem.MolFromSmiles(self.building_blocks[bb])
                            for bb in row[1:4]
                        ],
                        product_type="F",
                    )
                )
            except (ValueError, RuntimeError) as e:
                with self.subTest(vl_id=row[0], msg="Exception raised"):
                    self.fail(f"Error for vl_id {row[0]}: {e}")

            with self.subTest(vl_id=row[0]):
                self.assertEqual(
                    row[4],
                    generated_product,
                )

    @unittest.skipIf(not DB_EXISTS, "Missing database")
    def test_gives_correct_product_G(self):
        """
        Test that the reaction generator gives the correct product for a given set of starting materials
        on a large random subset of the VL.
        """

        res = self.con.con.execute(
            "SELECT id, initiator_long, monomer_long, terminator_long, SMILES FROM virtuallibrary WHERE type = 'G' AND initiator_long != '4-Pyrazole002';",
        ).fetchall()  # 4-Pyrazole002 is the bullshit initiator with "urea-KAT"

        for row in random.sample(res, 10000):  # we test 10000 random VL members
            try:
                generated_product = Chem.MolToSmiles(
                    self.rxn_generator.generate_product(
                        [
                            Chem.MolFromSmiles(self.building_blocks[bb])
                            for bb in row[1:4]
                        ],
                        product_type="G",
                    )
                )
            except (ValueError, RuntimeError) as e:
                with self.subTest(vl_id=row[0], msg="Exception raised"):
                    self.fail(f"Error for vl_id {row[0]}: {e}")

            with self.subTest(vl_id=row[0]):
                self.assertEqual(
                    row[4],
                    generated_product,
                )

    @unittest.skipIf(not DB_EXISTS, "Missing database")
    def test_gives_correct_product_H(self):
        """
        Test that the reaction generator gives the correct product for a given set of starting materials
        on a large random subset of the VL.
        """

        res = self.con.con.execute(
            "SELECT id, initiator_long, monomer_long, terminator_long, SMILES FROM virtuallibrary WHERE type = 'H' AND initiator_long != '4-Pyrazole002';",
        ).fetchall()  # 4-Pyrazole002 is the bullshit initiator with "urea-KAT"

        for row in random.sample(res, 10000):  # we test 10000 random VL members
            try:
                mol = self.rxn_generator.generate_product(
                    [Chem.MolFromSmiles(self.building_blocks[bb]) for bb in row[1:4]],
                    product_type="H",
                )
                if mol:
                    generated_product = Chem.MolToSmiles(mol)
                else:
                    # skip if the VL also has no product (None or "O" as placeholder)
                    if row[4] in ["O", None]:
                        continue
            except (ValueError, RuntimeError) as e:
                with self.subTest(vl_id=row[0], msg="Exception raised"):
                    self.fail(f"Error for vl_id {row[0]}: {e}")

            with self.subTest(vl_id=row[0], i=row[1], m=row[2], t=row[3]):
                self.assertEqual(
                    Chem.MolToSmiles(Chem.MolFromSmiles(row[4])),
                    generated_product,
                )

    def test_gives_correct_reaction_smiles(self):
        """
        Test that the reaction generator gives the correct reactionSMILES for a given product.
        """
        self.maxDiff = None
        # here we use a few manually curated samples to check if the result is correct
        products = [
            "Cc1ccccc1C(=O)NC1(Cc2nnc(-c3ccc4c(c3)OCO4)s2)CSC1",
            "CC(C)(C)OC(=O)n1c(-c2nnc(C3(CNC(=O)CCCCNC(=O)OCc4ccccc4)CCC3)s2)cc2ccccc21",
            "COc1cccc(C(=O)NCC2(c3nnc(-c4cn[nH]c4)s3)CCC2)n1",
        ]
        reaction_smiles = [
            "F[B-](F)(F)[C:1](=[O:2])[c:11]1[cH:12][cH:14][cH:17][cH:16][c:13]1[CH3:15].O=C1OC2(CCCCC2)O[C:6]12O[NH:3][C:4]1([CH2:5]2)[CH2:18][S:20][CH2:19]1.[C:7](=[S:8])([NH:9][NH2:10])[c:21]1[cH:22][cH:24][c:26]2[c:25]([cH:23]1)[O:27][CH2:29][O:28]2>>[C:1](=[O:2])([NH:3][C:4]1([CH2:5][c:6]2[s:8][c:7](-[c:21]3[cH:22][cH:24][c:26]4[c:25]([cH:23]3)[O:27][CH2:29][O:28]4)[n:9][n:10]2)[CH2:18][S:20][CH2:19]1)[c:11]1[cH:12][cH:14][cH:17][cH:16][c:13]1[CH3:15]",
            "F[B-](F)(F)[C:1](=[O:2])[CH2:11][CH2:12][CH2:13][CH2:14][NH:15][C:16](=[O:17])[O:18][CH2:19][c:20]1[cH:21][cH:23][cH:25][cH:24][cH:22]1.O=C1OC2(CCCCC2)O[C:6]12O[NH:3][CH2:4][C:5]21[CH2:26][CH2:28][CH2:27]1.[C:7](=[S:8])([NH:9][NH2:10])[c:29]1[n:30]([C:32]([O:35][C:39]([CH3:42])([CH3:43])[CH3:44])=[O:36])[c:33]2[c:34]([cH:31]1)[cH:38][cH:41][cH:40][cH:37]2>>[C:1](=[O:2])([NH:3][CH2:4][C:5]1([c:6]2[s:8][c:7](-[c:29]3[n:30]([C:32]([O:35][C:39]([CH3:42])([CH3:43])[CH3:44])=[O:36])[c:33]4[c:34]([cH:31]3)[cH:38][cH:41][cH:40][cH:37]4)[n:9][n:10]2)[CH2:26][CH2:28][CH2:27]1)[CH2:11][CH2:12][CH2:13][CH2:14][NH:15][C:16](=[O:17])[O:18][CH2:19][c:20]1[cH:21][cH:23][cH:25][cH:24][cH:22]1",
            "F[B-](F)(F)[C:1](=[O:2])[c:11]1[cH:12][cH:14][cH:16][c:15]([O:17][CH3:18])[n:13]1.O=C1OC2(CCCCC2)O[C:6]12O[NH:3][CH2:4][C:5]21[CH2:19][CH2:21][CH2:20]1.[C:7](=[S:8])([NH:9][NH2:10])[c:22]1[cH:23][n:25][nH:26][cH:24]1>>[C:1](=[O:2])([NH:3][CH2:4][C:5]1([c:6]2[s:8][c:7](-[c:22]3[cH:23][n:25][nH:26][cH:24]3)[n:9][n:10]2)[CH2:19][CH2:21][CH2:20]1)[c:11]1[cH:12][cH:14][cH:16][c:15]([O:17][CH3:18])[n:13]1",
        ]

        for prod, rsmi in zip(
            products, reaction_smiles
        ):  # we test 10000 random VL members
            try:
                generated_reaction = self.rxn_generator.get_reaction_smiles(prod)

            except (ValueError, RuntimeError) as e:
                with self.subTest(product=prod, msg="Exception raised"):
                    self.fail(f"Exception raised: {e}")

            with self.subTest(product=prod):
                self.assertEqual(
                    rsmi,
                    generated_reaction,
                )


if __name__ == "__main__":
    unittest.main()
