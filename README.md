# semiinf
## How to use semiinf
1. Edit make.inc and compile semiinf.x
	* make clean; make
2. (If bulk-slab matching is required) Edit input parameters in input_file.txt and run.
	* python3 util_match/match_main.py input_file.txt
3. Edit input file and do semiinfinite calculation
	* PATH/TO/SEMIINF/semiinf.x input_file.in > output_file.out

## How to use the wannier\_tb.py utility
wannier\_tb.py is a python script to read, process, create, and write tight-binding models in a wannier90-readable format (i.e. seedname\_hr.dat format).

In the semiinf program, wannier_tb.py is used to create \_hr.dat files for simple tight-binding models which can then be processed in semiinf.x to calculate surface and bulk DOS.
