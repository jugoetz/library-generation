# Library Generation Tools for Synthetic Fermentation

<a href="https://github.com/jugoetz/library-generation/blob/master/LICENSE"><img alt="License: MIT" src="https://black.readthedocs.io/en/stable/_static/license.svg"></a>
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>

A collection of tools needed to deal with HTE compound libraries.
Specialized on Synthetic Fermentation.
We used this code to plan, execute, and evaluate experiments in the 50k project (conducting and analyzing 50,000 chemical reactions).

For code used to conduct machine learning experiments, see [this repository](https://github.com/jugoetz/synferm-predictions).

## Installation

```bash
conda env create -f environment.yml
```

### Development
We use [nbstripout](https://pypi.org/project/nbstripout/) to remove output from notebooks before committing to the repository.
Install with:
```bash
conda install -c conda-forge nbstripout # or pip install nbstripout
nbstripout --install  # configures git filters and attributes for this repo
```
We use [pre-commit](https://pre-commit.com/) hooks to ensure standardized code formatting.
Install with:
```bash
pre-commit install
```

## Usage

General usage is as follows:

```bash
python -m src.path.to.submodule <args>
```

The configuration file is located in `src/config.yaml`.
Some modules accept command-line arguments, which, if given, override the configuration file settings.

### Submodules
The code is split over multiple submodules, each covering a different aspect of the library generation process.
Example workflows are given in the [Usage scenarios](#usage-scenarios) section below.

#### `analysis`

This submodule contains tools to evaluate experiments.
This means processing LC-MS data and instrument logs.

- `calculatelcmsyields`: Convert LC-MS results to ratios against internal standard.
- `extracterrors`: Extract errors from MoBiAS output, NEXUS logs.
- `extractmobiasresults`: Extract LC-MS results from MoBiAS output.
- `heatmapfromdb`: Generate yield heatmaps from processed LC-MS data.

#### `experiment_planning`

This submodule contains tools to set up individual experiments.

- `generatelcmssubmission`: Generate the submission file for MoBiAS LC-MS.
- `weighincalculator`: Produce a weigh-in sheet for a given experiment to prepare the stock solutions.
- `writereactiondb`: Populate the reaction database with a given experiment (reactants, location on the plate, ...).

#### `library_design`

This submodule contains tools to set up the library synthesis agenda.
This entails collecting building block data saved in an external inventory system,
cleaning and pruning the building block list,
enumerating the virtual library,
creating the overarching synthesis plan,
and individual plate layouts.

- `cheminventorycleanup`: Prune the building blocks stored in [ChemInventory](https://www.cheminventory.net/)
    to remove duplicates, scarce, and impure material.
    This requires the ChemInventory export file in `data/inputs/IventoryExport.xlsx`.
- `enumeratelibrary`: Take the pickled DataFrame from `inventorytosdf.py` and enumerate the corresponding virtual
  library. Save to SDF.
- `generate_source_transfer_files`: Generate transfer files for OT-2 to prepare source plates that can be used with an Echo.
- `generatelibraryplan`: Take the building block list received from `cheminventorycleanup.py` and
draw random combinations for every experiment to give a synthesis plan.
- `generateplatelayout`: Take the synthesis plan given by `generatelibraryplan.py` and generate the
  reaction plate layouts for every experiment.
- `inventorytosdf`: Take the building block list received from `cheminventorycleanup.py` and amend it with
  additional information (e.g. MW, weigh-in). Gives and SDF of building blocks and pickled DataFrame with the building block Mols.
- `layoutsourceplate`: From the synthesis plan generated in  `generatelibraryplan.py`, generate the source plate layouts for the 50k experiments.
- `layoutsourceplate_validation`: (this is a specialized script only for the 1D validation plates)
  From the synthesis plan for 1D validation plates, generate the source plate layouts.
- `reaction_generator`: Implements a class `SFReactionGenerator` that can be used to generate Synthetic Fermentation
    reactants from products, products from reactants, and atom-mapped reactionSMILES.
- `sdftoproperties`: Take the VL (as .sdf.gz files) and extract / calculate identifiers, molecular weights and
  formulae of the library constituents.

#### `utils`
Contains a range of utilities that we don't describe in detail here.

### Usage scenarios
#### When additional building blocks are added to the library

Two option:
1. If you want to remake the entire database (probably you don't), put the ChemInventory export file into <root>/data/inputs/ and run library_design scripts everything.
2. If you just want to add a few buildingblocks to the existing database, use `SynfermDatabaseConnection.add_building_block()` from `src.util.db_utils` to add the building blocks to the database.

In any case, you will next want to add new entries to the `virtuallibrary` table to reflect the updated building blocks.
To do so, run the Jupyter notebook `notebooks/virtuallibrary/add_new_products_to_vl.ipynb`.
This will determine combinations missing from the virtual library and add them.

#### When setting up a new experiment

1. Generate target plate layouts. For this step, there is no general procedure, but some tools that can be helpful:
  - `generateplatelayout.py` can be used to generate the target plate layouts from a synthesis plan (50k project).
  - `layoutsourceplate.py` can be used to generate the source plate layouts from the target plate layouts (50k project).
  - If source plate layout and transfer files are available, the `labware.transfers.Transfer` class has a
      `simulate()` method that can be used to generate the target plate layouts.
  - The Jupyter notebooks `notebooks/experiment_design/irregular_Echo_transfer_files_validation.ipynb` and
      `notebooks/experiment_design/validation_plates.ipynb` were used to generate non-standard source and transfer files.
2. If any building blocks are used that are not already in the database, follow the instructions above to populate the
    building_blocks and virtuallibrary tables with the new building blocks.
3. Use `src.experiment_planning.writereactiondb` to populate the reaction database with the new experiment(s).
    You only need the target plate layouts (step 1) and the populated virtuallibrary table (step 2).
4. Use `src.experiment_planning.weighincalculator` to generate the weigh-in sheet (not tested for general use, should work if volumes are adjusted in config file)
5. Use `src.experiment_planning.generatelcmssubmission` to generate the MoBiAS submission file.
    You only need the target plate layouts (step 1) and the populated experiments table (step 3).

#### When evaluating an experiment
At this point you should have received the MoBiAS output files and the NEXUS log files.
1. Place MoBiAS output files in the plate directory (e.g. `plates/JG404`).
2. Place NEXUS log files in the experiment directory (e.g. `plates/exp29`). Use the structure of `data/util/transfer_files_folder_template` as a guide.
3. Run `src.analysis.extractmobiasresults` to extract the raw LC-MS data from the MoBiAS output files
4. Run `src.analysis.extracterrors` to extract the errors from the MoBiAS output files and the NEXUS log files
5. Run `src.analysis.calculatelcmsyields` to calculate the yields from the raw LC-MS data
6. Run `src.analysis.heatmapfromdb` to generate the yield heatmaps from the processed LC-MS data.


#### When new products are added in the enumeration (legacy scenario, superseded by enumeration through `SFReactionGenerator` and `add_new_products_to_vl.ipynb`)
Change `enumerate_reaction.py` and rerun:

- `enumerate_reaction.py`
- `sdftoproperties.py`
- `generatelcmssubmission.py`

#### When the synthesis plan has to be changed (e.g. compound cannot be located)

- Apply changes to `synthesis_plan.json` manually.
  You should also do the changes in `synthesis_plan.csv`, but this will have no downstream influence.
- Run `generateplatelayout.py` to generate the new target plate layouts (this will operate in the `new` directory).
- Run `layoutsourceplate.py` to generate the source plate layouts (this will operate in the `new` directory).
- Manually copy appropriate directories from `new` to `plates`.

You are good to go.
Depending on previous work you might additionally have to delete now obsolete entries in the experiment database table,
re-rerun `writereactiondb.py` and `weighin_calculator.py`, and delete entries in `plate_list.csv`.
