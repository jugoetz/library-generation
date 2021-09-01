import re
from pathlib import Path

# directories and files
DATA_DIR = Path(
    '/Users/julian/PycharmProjects/library-generation/data').resolve()  # TODO find a better way to write this path as relative
INPUT_DIR = DATA_DIR / 'inputs'
BB_DIR = DATA_DIR / 'buildingblocks'
LIB_STATIC_DIR = DATA_DIR / 'library_static'
LIB_SDF_DIR = DATA_DIR / 'library_sdf'
LIB_INFO_DIR = DATA_DIR / 'library_info'
LOG_DIR = DATA_DIR / 'logs'
PLATES_DIR = DATA_DIR / 'plates'
DB_PATH = DATA_DIR / 'db' / '50k_project.db'
COMPOUND_MAPPING_PATH = DATA_DIR / 'buildingblocks' / 'compound_mapping.txt'
# COMPOUND_MAPPING_PATH = EXP_DIR / 'identity.txt'
# PLATE_REGEX = re.compile('test_JG([0-9]+).csv')
PLATE_REGEX = re.compile('plate_layout_plate([0-9]+).csv')

###################
EXP_NR = 'JG258'  # this will change frequently
EXP_DIR = PLATES_DIR / EXP_NR
