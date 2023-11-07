import unittest
import random

from rdkit import Chem

from src.library_design.reaction_generator import SFReactionGenerator
from src.util.db_utils import SynFermDatabaseConnection
from src.util.rdkit_util import desalt_building_block, remove_monomer_pg_chirality


class TestSFReactionGenerator(unittest.TestCase):
    def setUp(self) -> None:
        self.rxn_generator = SFReactionGenerator()
        self.con = SynFermDatabaseConnection()
        self.spiro_c_pattern = Chem.MolFromSmarts(
            "[$([CR2](O1)(ONC2)(C2)C(=O)OC1)]"
        )  # hits a spiro carbon atom between an isoxazolidine and a 5-membered lactone-ether ring
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
        }  # we only do this sanitization once per building block - much faster when we run the entire vl
        self.building_blocks[
            "TerTH010"
        ] = "NNC(=S)C=Cc1ccccc1"  # known error for double bond in VL

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

    def test_gives_correct_product(self):
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
            "F[B-](F)(F)[C:2]([c:1]1[cH:13][cH:15][cH:18][cH:17][c:14]1[CH3:16])=[O:3].O=C1OC2(CCCCC2)O[C:7]12O[NH:4][C:5]1([CH2:6]2)[CH2:19][S:21][CH2:20]1.[c:8]1([C:9](=[S:10])[NH:11][NH2:12])[cH:22][cH:24][c:26]2[c:25]([cH:23]1)[O:27][CH2:29][O:28]2>>[c:1]1([C:2](=[O:3])[NH:4][C:5]2([CH2:6][c:7]3[s:10][c:9](-[c:8]4[cH:22][cH:24][c:26]5[c:25]([cH:23]4)[O:27][CH2:29][O:28]5)[n:11][n:12]3)[CH2:19][S:21][CH2:20]2)[cH:13][cH:15][cH:18][cH:17][c:14]1[CH3:16]",
            "F[B-](F)(F)[C:2]([CH2:1][CH2:13][CH2:14][CH2:15][NH:16][C:17](=[O:18])[O:19][CH2:20][c:21]1[cH:22][cH:24][cH:26][cH:25][cH:23]1)=[O:3].O=C1OC2(CCCCC2)O[C:7]12O[NH:4][CH2:5][C:6]21[CH2:27][CH2:29][CH2:28]1.[c:8]1([C:9](=[S:10])[NH:11][NH2:12])[n:30]([C:32]([O:35][C:39]([CH3:42])([CH3:43])[CH3:44])=[O:36])[c:33]2[c:34]([cH:31]1)[cH:38][cH:41][cH:40][cH:37]2>>[CH2:1]([C:2](=[O:3])[NH:4][CH2:5][C:6]1([c:7]2[s:10][c:9](-[c:8]3[n:30]([C:32]([O:35][C:39]([CH3:42])([CH3:43])[CH3:44])=[O:36])[c:33]4[c:34]([cH:31]3)[cH:38][cH:41][cH:40][cH:37]4)[n:11][n:12]2)[CH2:27][CH2:29][CH2:28]1)[CH2:13][CH2:14][CH2:15][NH:16][C:17](=[O:18])[O:19][CH2:20][c:21]1[cH:22][cH:24][cH:26][cH:25][cH:23]1",
            "F[B-](F)(F)[C:2]([c:1]1[cH:13][cH:15][cH:17][c:16]([O:18][CH3:19])[n:14]1)=[O:3].O=C1OC2(CCCCC2)O[C:7]12O[NH:4][CH2:5][C:6]21[CH2:20][CH2:22][CH2:21]1.[c:8]1([C:9](=[S:10])[NH:11][NH2:12])[cH:23][n:25][nH:26][cH:24]1>>[c:1]1([C:2](=[O:3])[NH:4][CH2:5][C:6]2([c:7]3[s:10][c:9](-[c:8]4[cH:23][n:25][nH:26][cH:24]4)[n:11][n:12]3)[CH2:20][CH2:22][CH2:21]2)[cH:13][cH:15][cH:17][c:16]([O:18][CH3:19])[n:14]1",
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
