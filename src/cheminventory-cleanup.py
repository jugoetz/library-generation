
"""

This script performs cleanup operations on the building block libraries obtained from ChemInventory.

Inputs:
    - InventoryExport.xlsx: ChemInventory output Excel sheet
    - exceptions.txt: Space-delimited list of building blocks that should be ignored

The following entries are removed:
    - Chemicals not relevant to the project (e.g. SnAP)
    - empty bottles (N/A or 0 mass)
    - bottles below mass threshold of 10 mg
    - bottles/compounds specified in external list (specialty building blocks like dyes)
    - bottles carrying a tag "impure" in their comments
    - duplicates
The script further removes unneeded information (e.g. Supplier, CAS, References, ...).

Outputs:
    - inventory_containers.csv: List of individual containers (there can be multiple containers per compound)
    - inventory_compounds.csv: List of compounds with name, SMILES, MW, and category (I/M/T)
    - removed_small_amount.csv, removed_dimers.csv, removed_impure.csv: List of compounds removed at the respective step of
        the cleanup process

"""

import pandas as pd
from pathlib import Path

DATA_DIR = Path('..', 'data').resolve()

# import data
file = DATA_DIR / 'inputs' / 'InventoryExport.xlsx'
df = pd.read_excel(file, engine='openpyxl')

# import manual exceptions
with open(DATA_DIR / 'manual_settings' / 'exceptions.txt', 'r') as f:
    exceptions = [line.split(',')[0].strip('\n')
                  for line in f.readlines()
                  if line.startswith('#') is False and line != '\n'
                  ]
# ensure exceptions is not empty (otherwise it will delete entire df downstream)
if not exceptions:
    exceptions = ['PLACEHOLDER: No exceptions specified']


# drop irrelevant information
df.drop(columns=['Substance CAS',
                 'Supplier',
                 'Date Acquired',
                 'Molecular Formula',
                 'Molecular Weight',
                 'References',
                 'References.1',
                 'KAT',
                 'Class',
                 ],
        inplace=True
        )

# drop chemicals not needed (everything except I, M, T)
df.drop(df.loc[~df['Location'].str.contains('KATs|Monomers|Aminobenzenethiol|Thiohydrazide', regex=True), :].index, inplace=True)


# check if MW is not N/A
if df['MW [g/mol]'].isnull().sum() != 0:
    print('Warning: MW is not set for some bottles')


df['MW [g/mol]'] = pd.to_numeric(df['MW [g/mol]'])
df['Container Size'] = df['Container Size'].str.strip('<> ')
df['Container Size'] = pd.to_numeric(df['Container Size'])

# remove all bottles with mass N/A
df.dropna(axis=0, subset=['Container Size'], inplace=True)

# remove all bottles with mass <= 10 mg
df.loc[(df['Container Size'] <= 10) & (df['Unit'] == 'mg')].loc[:,'Container Name']\
    .to_csv(DATA_DIR / 'logs' / 'removed_small_amount.csv', index=False)
df.drop(df.loc[(df['Container Size'] <= 10) & (df['Unit'] == 'mg')].index, axis=0, inplace=True)

# remove bottles with "Dimer" as part of their name
df.loc[df['Container Name'].str.contains('Dimer', regex=False, case=False)].loc[:, 'Container Name']\
    .to_csv(DATA_DIR / 'logs' / 'removed_dimers.csv', index=False)
df.drop(df.loc[df['Container Name'].str.contains('Dimer', regex=False, case=False), :].index, axis=0, inplace=True)

# remove bottles that are specified in the exceptions list
df.drop(df.loc[df['Container Name'].str.contains('|'.join(exceptions), regex=True)].index, axis=0, inplace=True)

# remove bottles with "impure" in the comments
df.loc[df['Comments'].str.contains('impure', regex=False, na=False, case=False)].loc[:, 'Container Name']\
    .to_csv(DATA_DIR / 'logs' / 'removed_impure.csv', index=False)
df.drop(df.loc[df['Comments'].str.contains('impure', regex=False, na=False, case=False)].index, axis=0, inplace=True)

# introduce new column 'Compound Name' without trailing supplier initials or numbers
df['Compound Name'] = df['Container Name'].str.split('_', expand=True)[0]

# introduce new column 'Category' to distinguish I, M, T
df['Category'] = 'I'
df.loc[df['Location'].str.contains('Monomers', regex=False), 'Category'] = 'M'
df.loc[df['Location'].str.contains('Terminator', regex=False), 'Category'] = 'T'

# group by Compound Name
df_compounds = df.drop_duplicates(subset=['Compound Name'])
df_compounds = df_compounds[['Compound Name', 'SMILES', 'MW [g/mol]', 'Category']]

# report data
freq = df_compounds["Category"].value_counts()
print(f'The Virtual Library contains {freq["I"]} Initiators, {freq["M"]} Monomers, {freq["T"]} Terminators.\n'
      f'This results in a total of {freq["I"] * freq["M"] * freq["T"]:,} possible products.')

# save results
df.to_csv(DATA_DIR / 'outputs' / 'inventory_containers.csv', index=False)
df_compounds.to_csv(DATA_DIR / 'outputs' / 'inventory_compounds.csv', index=False)
