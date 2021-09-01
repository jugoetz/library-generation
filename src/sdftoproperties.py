"""
Take the gzipped sdf files that represent the VL and extract / calculate the name, molecular formula,
and monoisotopic mass of the library constituents from it. Dump the resulting dictionary to JSON for re-use
by other scripts.

The intent is to run this script on a powerful machine and make the extracted data available to downstream
scripts on weaker computers.

Needs to be rerun after changes to the static VL.
"""

import os
import gzip
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcMolFormula
import json
from config import *

"""GLOBALS"""
OUTPUT_JSON = LIB_STATIC_DIR / 'static_mol_prop_dict.json.gz'


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
    mol_prop_dict = {} # will hold the data
    files_to_process = []  # a helper to follow progress

    for path, _, files in os.walk(LIB_STATIC_DIR):
        # find which files we are interested in
        for f in files:
            if f.startswith('product_') and f.endswith('.sdf.gz'):
                files_to_process.append(f)
        # extract data from those files
        for i, f in enumerate(files_to_process):
            print(f'Now converting file {f}... [{i + 1}/{len(files_to_process)}]')
            letter = f.split('_')[1].split('.')[0]
            supplier = import_mol(Path(path, f))
            mol_props = {}
            for mol in supplier:
                if mol is not None:
                    mol_props[mol.GetProp('_Name')] = [CalcExactMolWt(mol), CalcMolFormula(mol)]
            mol_prop_dict[letter] = mol_props
            print(f'>> Done converting file {f}')

    # save
    print('All files converted. Saving to JSON...')
    with gzip.open(OUTPUT_JSON, 'wt', encoding='ascii') as zipfile:
        json.dump(mol_prop_dict, zipfile)
    print('### DONE! ###')

