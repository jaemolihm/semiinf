# semiinf examples
Each folder contains example input file for the `semiinf.x` program.
The data files (`\_hr.dat`, `\_spnr.dat`) are included in the `examplexx/reference\_out/` folder.
Also, one can open the ipython notebook `examplexx/examplexx.ipynb` and follow the pre- and post-processing scripts.

The tight-binding models used in the examples are generated either using the `wannier\_tb.py` utility, or using the wannier90.x program after DFT calculations.

## List of examples
- example01: s orbitals on a cubic lattice. Bulk and surface band structures and DOS.
 - Tight-binding model is generated using `util/wannier\_tb.py`
