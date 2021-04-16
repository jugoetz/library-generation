"""
Formerly first part of LCMS submission

Take the gzipped sdf files and construct the Name, molecular formula, and monoisotopic mass of the molecules from it
"""

import os
from pathlib import Path
import gzip
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcMolFormula
import pickle as pkl


"""GLOBALS"""
DATA_DIR = Path('..', 'data').resolve()
OUTPUT_DIR = DATA_DIR / 'outputs'
SDF_DIR = DATA_DIR / 'library_static'
OUTPUT_PICKLE = SDF_DIR / 'static_mol_prop_dict.pkl'
verbose = True


def import_mol(file):
    """
    Import product mol from gzipped SDF. Return a generator that yields rdkit mol objects
    :param file: str or Path-like
    :return: generator
    """
    if str(file).endswith('.sdf.gz'):
        supply = ForwardSDMolSupplier(gzip.open(file))
    else:
        raise ValueError(f'Invalid file "{file}". User must supply .sdf.gz file.')
    return supply


########################################################################################################################
#                                                        MAIN                                                          #
########################################################################################################################

if __name__ == '__main__':
    """Generate mol suppliers to import from SDF"""
    mol_prop_dict = {}
    for path, _, files in os.walk(SDF_DIR):
        for f in files:
            if f.startswith('product_') and f.endswith('.sdf.gz'):
                print(f'Now converting file {f}...')
                letter = f.split('_')[1].split('.')[0]
                supplier = import_mol(Path(path, f))
                mol_props = {}
                for mol in supplier:
                    if mol is not None:
                        mol_props[mol.GetProp('_Name')] = [CalcExactMolWt(mol), CalcMolFormula(mol)]
                mol_prop_dict[letter] = mol_props
                print(f'>> Done converting file {f}')
    print('All files converted. Saving to pickle...')
    with open(OUTPUT_PICKLE, 'wb') as file:
        pkl.dump(mol_prop_dict, file)
    print('### DONE! ###')
