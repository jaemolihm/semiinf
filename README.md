# semiinf
## How to use semiinf
1. Edit make.inc and compile semiinf.x
	* make clean; make
2. (If bulk-slab matching is required) Edit input parameters in input_file.txt and run.
	* python3 util_match/match_main.py input_file.txt
3. Edit input file and do semiinfinite calculation
	* PATH/TO/SEMIINF/semiinf.x input_file.in > output_file.out
