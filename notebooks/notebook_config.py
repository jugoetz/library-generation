import re
from pathlib import Path

# directories and files
DATA_DIR = Path('..', 'data').resolve()
INPUT_DIR = DATA_DIR / 'inputs'
SDF_DIR = DATA_DIR / 'library_static'
DB_PATH = DATA_DIR / 'db' / '50k_project.db'
DB_STATIC_DIR = DATA_DIR / 'db' / 'static'
SUBMISSION_FORM_TEMPLATE = DATA_DIR / 'util' / 'BMIIyyyyyy-SampleTable_JGxxx.xls'
EXP_NR = 'JG216'  # this will change frequently
EXP_DIR = DATA_DIR / 'plates' / EXP_NR

# debugging
DEBUG = 1
verbose = True

# Library things
# COMPOUND_MAPPING_PATH = EXP_DIR / 'identity.txt'
COMPOUND_MAPPING_PATH = DATA_DIR / 'buildingblocks' / 'compound_mapping.txt'
# PLATE_REGEX = re.compile('test_JG([0-9]+).csv')
PLATE_REGEX = re.compile('plate_layout_plate([0-9]+).csv')

# Chemistry things
ADD_IS = 'y'  # Was internal standard added to the plates?
IS_FORMULA = "C20H21O4Cl"  # molecular formula of Fenofibrat
IS_MASS = 360.1128  # monoisotopic mass of Fenofibrat
MASS_OR_FORMULA = 'formula'  # ['mass'/'formula'] can output mass as a number or can give chemical formula
PLATE_SIZE = 384
