## Library Generation Tools

Tools needed to deal generate HTE compound libraries. Specialized on Synthetic Fermentation.

### Submodules

#### db_interconversion

- **calculate_lcms_yields.py**: Conversion from raw LCMS data to ratios against internal standard

#### db_models

ORM models for database. Currently not use din production

#### db_retrieval

Tools to retrieve data from db, e.g., for visualization, or experiment planning

- **generatelcmssubmission.py**: Generate the submission file for Mobias
- **heatmapfromdatabase.py**: Obtain a visual representation of results from LCMS peak areas

#### db_writing

Addition of raw and processed data to database

#### legacy_code

as the name suggests...

#### library_design

Tools for going from building blocks save in an external inventory system to an experiment library

- **cheminventory-cleanup.py**: Prune the building blocks stored in cheminventory to remove duplicates, scarce, and
  impure material

- **inventorytosdf.py**: Take the building block list received from cheminventory_cleanup.py and amend it with
  additional information (e.g. MW, weigh-in). Gives sdf of building blocks and serializes DataFrame of mols to pickle

- **enumerate-library.py**: Take the pickled DataFrame from inventorytosdf.py and enumerate the corresponding virtual
  library. Save to sdf.

- **generatelibraryplan.py**: Take the building block list received from cheminventory_cleanup.py and 
draw random combinations for every experiment to give a synthesis plan
  
- **generateplatelayout.py**: Take the synthesis plan given by generatelibraryplan.py and generate the 
target plates for every experiment.
  
- **sdf-to-properties.py**: Take the VL (as sdf.gz files) and extract / calculate names, molecular masses, and
formulae of the library molecules.
  
### Usage scenarios
#### When additional building blocks are added to the library

Put the ChemInventory export file into <root>/data/inputs/ and run library_design scripts everything

#### When new products are added in the enumeration

Change enumerate_reaction.py and rerun:

- enumerate_reaction.py
- sdf-to-properties.py
- generatelcmssubmission.py

#### When the synthesis plan has to be changed (e.g. compound cannot be located)

- Apply changes to `synthesis_plan.json` manually.
  You should also do the changes in `synthesis_plan.csv`, but this will have no downstream influence.
- Run `generateplatelayout.py` to generate the new target plate layouts (this will operate in the `new` directory).
- Run `layoutsourceplate.py` to generate the source plate layouts (this will operate in the `new` directory).
- Manually copy appropriate directories from `new` to `plates`.

You are good to go.
Depending on previous work you might additionally have to delete now obsolete entries in the experiment database table,
re-rerun `writereactiondb.py` and `weighin_calculator.py`, and delete entries in `plate_list.csv`.
