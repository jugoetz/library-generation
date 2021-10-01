from rdkit.Chem import Draw

from db_queries import MyDatabaseConnection

if __name__ == '__main__':
    lab_journal_number = 'JG263'
    well = 'H15'
    product = 'G'
    Draw.MolToImage(MyDatabaseConnection().get_product_mol_by_well(lab_journal_number, well, product),
                    size=(500, 500)).show()
