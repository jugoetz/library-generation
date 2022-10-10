import argparse

from rdkit.Chem import Draw, MolToSmiles

from src.util.db_utils import MyDatabaseConnection

defaults = {
    "lab_journal_nr": "JG290",
    "well": "A3",
    "product_type": "A",
    "output": "smiles",
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lab_journal_number",
        type=str,
        help="Lab journal number for the plate",
        default=defaults["lab_journal_nr"],
    )
    parser.add_argument(
        "--well",
        type=str,
        help="Well to show the product for",
        default=defaults["well"],
    )
    parser.add_argument(
        "--product_type",
        type=str,
        help="Product type to show",
        default=defaults["product_type"],
        choices=["A", "B", "C", "D", "E", "F", "G", "H"],
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output type",
        default=defaults["output"],
        choices=["smiles", "image"],
    )
    args = parser.parse_args()

    con = MyDatabaseConnection()

    if args.output == "image":
        Draw.MolToImage(
            con.get_product_mol_by_well(
                args.lab_journal_number, args.well, args.product_type
            ),
            size=(500, 500),
        ).show()
    elif args.output == "smiles":
        print(
            MolToSmiles(
                con.get_product_mol_by_well(
                    args.lab_journal_number, args.well, args.product_type
                )
            )
        )
