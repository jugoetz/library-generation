### How to use the library\_static directory

This directory holds sdf files generated on the MacPro to make them available for usage on
slower machines.

All files are gzipped for storage space reasons (use gzip -d [file] to decompress or use rdkit's capability to read gzipped sdfs directly).

Scripts should not write to this directory but may read from it.

### Library properties

Currently (12.03.2021), the library consists of 236,652 unique products from both TH and ABT reactions.
