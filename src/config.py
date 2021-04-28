import re
from pathlib import Path

# directories and files
DATA_DIR = Path('..', 'data').resolve()
OUTPUT_DIR = DATA_DIR / 'outputs'
INPUT_DIR = DATA_DIR / 'inputs'
SDF_DIR = DATA_DIR / 'library_static'
DB_PATH = DATA_DIR / 'db' / '50k_project.db'
EXP_NR = 'JG216'  # this will change frequently
EXP_DIR = OUTPUT_DIR / 'target_plates' / EXP_NR

# debugging
DEBUG = 1
verbose = True

# Library things
# COMPOUND_MAPPING = EXP_DIR / 'identity.txt'
COMPOUND_MAPPING = OUTPUT_DIR / 'compound_mapping.txt'
# PLATE_REGEX = re.compile('test_JG([0-9]+).csv')
PLATE_REGEX = re.compile('plate_layout_plate([0-9]+).csv')

# Chemistry things
ADD_IS = 'y'  # Was internal standard added to the plates?
IS_FORMULA = "C20H21O4Cl"  # molecular formula of Fenofibrat
IS_MASS = 360.1128  # monoisotopic mass of Fenofibrat
MASS_OR_FORMULA = 'formula'  # ['mass'/'formula'] can output mass as a number or can give chemical formula
PLATE_SIZE = 384
