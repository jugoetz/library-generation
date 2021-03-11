## Library Generation Tools

Tools needed to deal generate HTE compound libraries. Specialized on Synthetic Fermentation.

### Constituents
- **cheminventory-cleanup.py**: Prune the building blocks stored in cheminventory to remove duplicates,
scarce, and impure material
  
- **generatelibraryplan.py**: Take the building block list received from cheminventory_cleanup.ipynb and 
draw random combinations for every experiment to give a synthesis plan
  
- **generateplatelayout.py**: Take the synthesis plan given by generatelibraryplan.py and generate the 
target plates for every experiment. 

### Next steps

- from the target plate layouts, generate the analysis data for Mobias and the transfer files for Nexus