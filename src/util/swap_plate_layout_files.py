"""
This script is part of the error correction for the error discovered on 11.08.2022 related to the incorrect transfer
file. The error would cause plate Synthesis 3 to contain what was meant to go into Synthesis 1. Likewise, plate
Synthesis1 contains what was meant to go into Synthesis2 and plate Synthesis2 contains what was meant to go into
Synthesis3.

As part of the correction, this script will take an experiment directory, backup the plate layout files, then swap
them according to
    Synthesis1 -> Synthesis3
    Synthesis2 -> Synthesis1
    Synthesis3 -> Synthesis2
"""
import shutil

from src.definitions import PLATES_DIR

#################################
# set this directory before use #
exp_dir = PLATES_DIR / "exp18"  #
#################################

filenames = [f"plate_layout_plate{i + 1}.csv" for i in range(3)]
plate_layout_files = [exp_dir / f for f in filenames]
backup_dir = exp_dir / "legacy_2022-08-12"

try:
    backup_dir.mkdir()
except FileExistsError:
    raise FileExistsError(f"Backup directory {backup_dir} already exists. Please remove it and try again.")

for file in plate_layout_files:
    shutil.copy(file, backup_dir)

# rename the files
temp_file = plate_layout_files[0].rename(exp_dir / "temp.csv")
plate_layout_files[1].rename(exp_dir / "plate_layout_plate1.csv")
plate_layout_files[2].rename(exp_dir / "plate_layout_plate2.csv")
temp_file.rename(exp_dir / "plate_layout_plate3.csv")
