from pathlib import Path

ROOT_DIR = Path(__file__).parents[1]

DATA_DIR = ROOT_DIR / 'data'

INPUT_DIR = DATA_DIR / 'inputs'
DB_DIR = DATA_DIR / 'db'
BUILDING_BLOCKS_DIR = DATA_DIR / 'buildingblocks'
IMAGES_DIR = DATA_DIR / 'images'
LIB_INFO_DIR = DATA_DIR / 'library_info'
LIB_SDF_DIR = DATA_DIR / 'library_sdf'
LIB_STATIC_DIR = DATA_DIR / 'library_static'
LOG_DIR = DATA_DIR / 'logs'
MANUAL_SETTINGS_DIR = DATA_DIR / 'manual_settings'
PLATES_DIR = DATA_DIR / 'plates'
UTIL_DIR = DATA_DIR / 'util'

COMPOUND_MAPPING_PATH = BUILDING_BLOCKS_DIR / 'compound_mapping.txt'
# COMPOUND_MAPPING_PATH = EXP_DIR / 'identity.txt'   # this is sometimes useful for non-canonical experiments
CONF_PATH = ROOT_DIR / 'src' / 'config.yaml'
DB_PATH = DB_DIR / '50k_project.db'
PLATE_LIST_PATH = PLATES_DIR / 'plate_list.csv'
