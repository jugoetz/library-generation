"""

Script to amend the data from the inventory (after cleanup) with additional information

Inputs:
    - inventory_compounds.csv: List of library constituents with Compound Name, SMILES code, and category

Steps:
    - SMILES -> RdKit MOL -> image
    - MW -> weigh-in in mg for 100 µL of a 0.05 M solution (with rounding)
    - Resort columns for convenience

- Outputs:
    - all.sdf: Structures and info for all molecules
    - initiators.sdf, monomers.sdf, terminators.sdf: Structures and info divided by building blocks
    - Write an excel sheet in the style of inventory_compounds.csv with added picture and weigh-in
      (basically a excel dump of the df)
    - Pickle the df to reuse it easier

"""

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem.PropertyMol import PropertyMol
from pathlib import Path
import pandas as pd
import math
from config import *

"""Generate a DataFrame with all relevant information from inventory data"""
print('Importing inventory data...')
compounds = pd.read_csv(OUTPUT_DIR / 'inventory_compounds.csv')  # read df from inventory data
compounds['mol'] = compounds['SMILES'].apply(Chem.MolFromSmiles)  # generate rdkit mol objects
compounds['mol'] = compounds['mol'].apply(PropertyMol)
compounds.apply(lambda x: x['mol'].SetProp('_Name', x['Compound Name']), axis=1)  # add the name to mol for later saving to sdf
compounds['img'] = compounds['mol'].apply(Draw.MolToImage)  # generate molecule images
compounds['exact mass'] = compounds['mol'].apply(Descriptors.ExactMolWt)
compounds['MW_from_mol'] = compounds['mol'].apply(Descriptors.MolWt)
compounds['weigh-in [mg] / 100 µL'] = compounds['MW_from_mol']\
    .apply(lambda x: round(x * 1e-4 * 0.05 * 1000, 2))  # calculate the weigh-in mass from MW from structure

"""Control: Check if ChemInventory MW and MW from structure align"""
print('Verifying molecular weights...')
for i, data in compounds.iterrows():
    if not math.isclose(data['MW [g/mol]'], data['MW_from_mol'], abs_tol=0.1):
        print('WARNING: Molecular weight for the following compound is not in agreement with the value in ChemInventory')
        print(f'Cheminventory MW: {data["MW [g/mol]"]}')
        print(f'Calculate MW from mol: {data["MW_from_mol"]}')
        print(data)
compounds.drop(columns=['MW [g/mol]'], inplace=True)  # we won't use this. Calculation from structure is more reliable.

"""Generate outputs"""
# output to Excel
print('Generating Excel printout...')
compounds.drop(columns=['mol', 'img'], inplace=False).to_excel(OUTPUT_DIR / 'inventory_compounds_extended.xlsx')

# write the molecule images to files
print('Generating molecule images...')
for i, data in compounds.iterrows():
    with open(DATA_DIR / 'images' / ''.join([data.loc['Compound Name'], '.png']), 'wb') as file:
        data.loc['img'].save(file)

img = Draw.MolsToGridImage(compounds['mol'].tolist(),
                           molsPerRow=6,
                           subImgSize=(400, 400),
                           legends=compounds['Compound Name'].tolist(),
                           )
img_I = Draw.MolsToGridImage(compounds['mol'].loc[compounds['Category'] == 'I'].tolist(),
                           molsPerRow=4,
                           subImgSize=(600, 600),
                           legends=compounds['Compound Name'].loc[compounds['Category'] == 'I'].tolist(),
                           )
img_M = Draw.MolsToGridImage(compounds['mol'].loc[compounds['Category'] == 'M'].tolist(),
                           molsPerRow=4,
                           subImgSize=(600, 600),
                           legends=compounds['Compound Name'].loc[compounds['Category'] == 'M'].tolist(),
                           )
img_T = Draw.MolsToGridImage(compounds['mol'].loc[compounds['Category'] == 'T'].tolist(),
                           molsPerRow=4,
                           subImgSize=(600, 600),
                           legends=compounds['Compound Name'].loc[compounds['Category'] == 'T'].tolist(),
                           )

# write a big overview
with open(DATA_DIR / 'images' / '_overview.png', 'wb') as file:
    img.save(file)
with open(DATA_DIR / 'images' / '_I.png', 'wb') as file:
    img_I.save(file)
with open(DATA_DIR / 'images' / '_M.png', 'wb') as file:
    img_M.save(file)
with open(DATA_DIR / 'images' / '_T.png', 'wb') as file:
    img_T.save(file)

# output to SDF
print('Generating SDFs...')
with open(OUTPUT_DIR / 'sdf' / 'initiators.sdf', 'w') as file_i, \
        open(OUTPUT_DIR / 'sdf' / 'monomers.sdf', 'w') as file_m, \
        open(OUTPUT_DIR / 'sdf' / 'terminators.sdf', 'w') as file_t:
    writer_i = Chem.SDWriter(file_i)
    writer_m = Chem.SDWriter(file_m)
    writer_t = Chem.SDWriter(file_t)
    for i, data in compounds.iterrows():
        if data.loc['Category'] == 'I':
            writer_i.write(data.loc['mol'])
        if data.loc['Category'] == 'M':
            writer_m.write(data.loc['mol'])
        if data.loc['Category'] == 'T':
            writer_t.write(data.loc['mol'])
    """It is necessary to close the SDWriters manually to prevent an exception during garbage collection"""
    writer_i.close()
    writer_m.close()
    writer_t.close()

# dump df
print('Generating pickle dump...')
compounds.to_pickle(OUTPUT_DIR / 'library_constituents_dataframe.pkl')

print('\nFinished. Exiting...')
