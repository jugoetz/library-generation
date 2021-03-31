## Library Generation Tools

Tools needed to deal generate HTE compound libraries. Specialized on Synthetic Fermentation.

### Constituents
- **cheminventory-cleanup.py**: Prune the building blocks stored in cheminventory to remove duplicates,
scarce, and impure material
  
- **inventorytosdf.py**: Take the building block list received from cheminventory_cleanup.py and amend it with
additional information (e.g. MW, weigh-in). Gives sdf of building blocks and serializes DataFrame of mols to pickle
  
- **enumerate-library.py**: Take the pickled DataFrame from inventorytosdf.py and enumerate the corresponding
  virtual library. Save to sdf.
  
- **generatelibraryplan.py**: Take the building block list received from cheminventory_cleanup.py and 
draw random combinations for every experiment to give a synthesis plan
  
- **generateplatelayout.py**: Take the synthesis plan given by generatelibraryplan.py and generate the 
target plates for every experiment. 
  
- **generatelcmssubmission.py**: Take compound mapping, plate layout, and SDF files and generate the 
  submission file for Mobias

### Next steps

- from the target plate layouts, generate the transfer files for Nexus