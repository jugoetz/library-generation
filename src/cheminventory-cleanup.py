"""
This script performs cleanup operations on the building block libraries obtained from ChemInventory.

Inputs:
    - InventoryExport.xlsx: ChemInventory output Excel sheet
    - blacklist.txt: Space-delimited list of building blocks that should be ignored

Filters applied:
    - Chemicals not relevant to the project (e.g. SnAP)
    - empty bottles (N/A or 0 mass)
    - bottles below mass threshold of 10 mg
    - bottles/compounds specified in external list (specialty building blocks like dyes)
    - bottles carrying a tag "impure" in their comments
    - duplicates (how?)
The script further removes unneeded information (e.g. Supplier, CAS, References, ...). (not necessary)

Outputs:
    - inventory_containers.csv: List of individual containers (there can be multiple containers per compound)
    - inventory_compounds.csv: List of compounds with name, SMILES, MW, and category (I/M/T)
    - removed_small_amount.csv, removed_dimers.csv, removed_impure.csv: List of compounds removed at the respective step of
        the cleanup process

"""

import pandas as pd
from config import *


def import_from_xlsx():
    file = INPUT_DIR / 'InventoryExport.xlsx'
    df = pd.read_excel(file, engine='openpyxl')

    return df


def import_blacklist(path):
    with open(path, 'r') as f:
        exceptions = [line.split(',')[0].strip('\n')
                      for line in f.readlines()
                      if line.startswith('#') is False and line != '\n'
                      ]
    # ensure exceptions is not empty (otherwise it will delete entire df downstream)
    if not exceptions:
        exceptions = ['PLACEHOLDER: No exceptions specified']
    return exceptions


def filter_relevant_columns(df):
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
    return


def filter_relevant_categories(df):
    df.drop(df.loc[~df['Location'].str.contains('KATs|Monomers|Aminobenzenethiol|Thiohydrazide', regex=True), :].index,
            inplace=True)
    return


def convert_dtypes(df):
    # check if MW is not N/A
    if df['MW [g/mol]'].isnull().sum() != 0:
        print('Warning: MW is not set for some bottles')
    df['MW [g/mol]'] = pd.to_numeric(df['MW [g/mol]'])
    df['Container Size'] = df['Container Size'].str.strip('<> ')
    df['Container Size'] = pd.to_numeric(df['Container Size'])


def filter_mass(df):
    df.dropna(axis=0, subset=['Container Size'], inplace=True)
    df.loc[(df['Container Size'] <= 10) & (df['Unit'] == 'mg')].loc[:, 'Container Name'] \
        .to_csv(DATA_DIR / 'logs' / 'removed_small_amount.csv', index=False)
    df.drop(df.loc[(df['Container Size'] <= 10) & (df['Unit'] == 'mg')].index, axis=0, inplace=True)


def filter_dimer(df):
    df.loc[df['Container Name'].str.contains('Dimer', regex=False, case=False)].loc[:, 'Container Name'] \
        .to_csv(DATA_DIR / 'logs' / 'removed_dimers.csv', index=False)
    df.drop(df.loc[df['Container Name'].str.contains('Dimer', regex=False, case=False), :].index, axis=0, inplace=True)


def filter_blacklist(df):
    # import manual blacklist
    blacklist = import_blacklist(DATA_DIR / 'manual_settings' / 'blacklist.txt')  # TODO consider removing global
    df.drop(df.loc[df['Container Name'].str.contains('|'.join(blacklist), regex=True)].index, axis=0, inplace=True)


def filter_purity(df):
    df.loc[df['Comments'].str.contains('impure', regex=False, na=False, case=False)].loc[:, 'Container Name'] \
        .to_csv(DATA_DIR / 'logs' / 'removed_impure.csv', index=False)
    df.drop(df.loc[df['Comments'].str.contains('impure', regex=False, na=False, case=False)].index, axis=0,
            inplace=True)


def apply_filters(df):
    """Filters must act on df inplace and not give a return value"""
    # drop irrelevant information
    filter_relevant_columns(df)

    # drop chemicals not needed (everything except I, M, T)
    filter_relevant_categories(df)

    # convert MW and Container Size to numeric types
    convert_dtypes(df)

    # remove all bottles with mass N/A and <= 10 mg
    filter_mass(df)

    # remove bottles with "Dimer" as part of their name
    filter_dimer(df)

    # remove blacklisted bottles
    filter_blacklist(df)

    # remove bottles with "impure" in the comments
    filter_purity(df)

    return


df = import_from_xlsx()
apply_filters(df)

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
df.to_csv(BB_DIR / 'inventory_containers.csv', index=False)
df_compounds.to_csv(BB_DIR / 'inventory_compounds.csv', index=False)
