from rdkit.Chem import Draw, MolToSmiles

from db_queries import MyDatabaseConnection

con = MyDatabaseConnection()

if __name__ == '__main__':
    lab_journal_number = 'JG276'
    well = 'A5'

    product = 'A'
    how = 'smiles'

    if how == 'drawing':
        Draw.MolToImage(con.get_product_mol_by_well(lab_journal_number, well, product),
                        size=(500, 500)).show()
    elif how == 'smiles':
        print(MolToSmiles(con.get_product_mol_by_well(lab_journal_number, well, product)))
