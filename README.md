# Library Generation Tools

A collection of tools needed to deal generate HTE compound libraries. Specialized on Synthetic Fermentation.

## Installation

```bash
conda env create -f environment.yml
```
Place the ChemInventory export file into `data/inputs/IventoryExport.xlsx`.

## Usage

General usage is as follows:

```bash
python -m src.path.to.submodule <args>
```

The configuration file is located in `src/config.yaml`.
Some modules accept command-line arguments, which, if given, override the configuration file settings.

### Submodules

#### analysis

Tools to evaluate experiments.

- **calculatelcmsyields**: Conversion from raw LCMS data to ratios against internal standard
- **extracterrors**: Extracts errors from MoBiAS output, NEXUS logs.
- **extractmobiasresults**: Extracts raw LCMS data from MoBiAS output.
- **heatmapfromdb**: Generates yield heatmaps from processed LCMS data.

#### experiment_planning

Tools to set up individual experiments.

- **generatelcmssubmission**: Generate the submission file for MoBiAS.

- **weighincalculator**: Produce a weigh-in sheet for a given experiment
- **writereactiondb**: Populate the reaction database with a given experiment.

#### legacy_code

as the name suggests...

#### library_design

Tools for going from building blocks save in an external inventory system to an experiment library

- **cheminventorycleanup**: Prune the building blocks stored in cheminventory to remove duplicates, scarce, and
  impure material.

- **enumeratelibrary**: Enumerate the virtual library from a set of building blocks.

- **generateexperiment**: This is an experimental script for experiment design that is not used in production.

- **inventorytosdf**: Take the building block list received from cheminventory_cleanup.py and amend it with
  additional information (e.g. MW, weigh-in). Gives sdf of building blocks and serializes DataFrame of mols to pickle

- **enumeratelibrary**: Take the pickled DataFrame from inventorytosdf.py and enumerate the corresponding virtual
  library. Save to sdf.

- **generatelibraryplan**: Take the building block list received from cheminventory_cleanup.py and
draw random combinations for every experiment to give a synthesis plan

- **generateplatelayout**: Take the synthesis plan given by generatelibraryplan.py and generate the
  target plates for every experiment.

- **sdftoproperties**: Take the VL (as sdf.gz files) and extract / calculate names, molecular masses, and
  formulae of the library molecules.

#### utils
- **db_utils**: Contains DB connection class which defines many popular queries.
- **number_of_synthesized_compounds**: Calculates the number of synthesized compounds from the reaction database.
- **show_product**: Show the product(s) in a given well as SMILES or image

### Usage scenarios
#### When additional building blocks are added to the library

Put the ChemInventory export file into <root>/data/inputs/ and run library_design scripts everything

#### When new products are added in the enumeration

Change enumerate_reaction.py and rerun:

- enumerate_reaction.py
- sdftoproperties.py
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
