import sys
from time import localtime, strftime
import numpy as np
import numpy.matlib
import scipy.linalg

from match_module_read_data import (tb_read_data, tb_sort_wcenter, tb_shift_fermi,
    tb_correct_data)
from match_module_ham_process import run_opt_diagonly
from match_module_ham_mismatch import (ham_mismatch_intra, ham_mismatch_inter,
                                       ham_maximal_hopping, write_hr_dat)


################
# input parameters #
material = 'GeTe'
path = './'
input_params = {}

input_params['seedname_b'] = 'GeTe.bulk'
input_params['seedname_s'] = 'GeTe.slab'
input_params['isspinor'] = True

input_params['iatm_add'] = 4
input_params['nbnd_s'] = 310
input_params['nw_s'] = 248
input_params['nklist_s'] = np.array([7,7,1])
input_params['nbnd_b'] = 60
input_params['nw_b'] = 48
input_params['nklist_b'] = np.array([7,7,6])
input_params['iw_head'] = 104
input_params['iw_tail'] = 152
input_params['norbital_list'] = [8]*6

input_params['alat_s'] = np.array([[ 4.177748, 0.000000,  0.000000],
                   [-2.088874, 3.618036,  0.000000],
                   [ 0.000000, 0.000000, 62.869262]]).transpose()
input_params['alat_b'] = np.array([[ 4.177748, 0.000000,  0.000000],
                   [-2.088874, 3.618036,  0.000000],
                   [ 0.000000, 0.000000, 10.478210]]).transpose()
################################
# read binary files and initialize them
################################
print('begin initialization', flush=True)
tbbulk = tb_read_data(material=material, path=path, isbulk=True, input_params=input_params)
tbslab = tb_read_data(material=material, path=path, isbulk=False, input_params=input_params)
tb_sort_wcenter(tbbulk, tbslab)
# tb_correct_data(tbbulk, tbslab)
tb_shift_fermi(tbbulk, tbslab)
print('end initialization', flush=True)
################################
# for profiling
# import cProfile
# cProfile.run('run_opt_diagonly(tbbulk, tbslab, max_iter=200)')

# iden
ham_intra_max_iden, ham_intra_sqsum_iden = ham_mismatch_intra(tbslab['hr'], tbbulk, tbslab)
ham_inter_max_iden, ham_inter_sqsum_iden = ham_mismatch_inter(tbslab['hr'], tbbulk, tbslab)

# # ham_diag
obj_list, is_converged = run_opt_diagonly(tbbulk, tbslab, max_iter=10000)
ham_intra_max_ham, ham_intra_sqsum_ham = ham_mismatch_intra(tbslab['hr_ham'], tbbulk, tbslab)
ham_inter_max_ham, ham_inter_sqsum_ham = ham_mismatch_inter(tbslab['hr_ham'], tbbulk, tbslab)

ham_intra, ham_inter = ham_maximal_hopping(tbbulk, tbslab)

## printing
print("begin output", flush=True)

# save output hamiltonian as text file
write_hr_dat(tbslab)

# write result summary
with open(path+'out_match.txt', 'w') as f:
    f.write("########################################################\n")
    if not is_converged: f.write("#####WARNING: convergence not reached during minimization#####\n\n")
    f.write("Calculation main_match.py at " 
        + strftime("%a, %d %b %Y %H:%M:%S +0000", localtime())+"\n")
    f.write(f"Material: {tbslab['material']}\n")
    f.write(f"Raw data path: {tbbulk['path']}\n")
    f.write(f"Is spinor calculation: {tbslab['isspinor']}\n")
    f.write(f"nbnd_b = {tbbulk['nbnd']}, nbnd_s = {tbslab['nbnd']}\n")
    f.write(f"nw_b = {tbbulk['nw']}, nw_s = {tbslab['nw']}\n")
    f.write(f"nk_b = {tbbulk['nk']}, nk_s = {tbslab['nk']}\n")
    f.write("norbital_list : " + ' '.join(str(e) for e in tbbulk['norbital_list']) + "\n")
    f.write(f"maximal intra-layer hopping (eV) = {ham_intra:.10f}\n")
    f.write(f"maximal inter-layer hopping (eV) = {ham_inter:.10f}\n")
    f.write("\n")
    f.write(f"hamiltonian mismatch (corrected) (eV^2) = {obj_list[-1]:.10f}\n")
    f.write("\n")
    f.write("===hamiltonian mismatch, intra-layer, maximum error (eV)===\n")
    f.write(f"ham_intra_max_iden = {ham_intra_max_iden:.10f}\n")
    f.write(f"ham_intra_max_ham  = {ham_intra_max_ham:.10f}\n")
    f.write("\n")
    f.write("===hamiltonian mismatch, intra-layer, sum of square (eV^2)===\n")
    f.write(f"ham_intra_sqsum_iden = {ham_intra_sqsum_iden:.10f}\n")
    f.write(f"ham_intra_sqsum_ham  = {ham_intra_sqsum_ham:.10f}\n")
    f.write("\n")
    f.write("===hamiltonian mismatch, inter-layer, maximum error (eV)===\n")
    f.write(f"ham_inter_max_iden = {ham_inter_max_iden:.10f}\n")
    f.write(f"ham_inter_max_ham  = {ham_inter_max_ham:.10f}\n")
    f.write("\n")
    f.write("===hamiltonian mismatch, inter-layer, sum of square (eV^2)===\n")
    f.write(f"ham_inter_sqsum_iden = {ham_inter_sqsum_iden:.10f}\n")
    f.write(f"ham_inter_sqsum_ham  = {ham_inter_sqsum_ham:.10f}\n")
    f.write("########################################################\n")

# output minimization result in gnuplot format
with open(path+'out_match_minimization.txt', 'w') as f:
    f.write("# iter ham_mismatch(eV^2)\n")
    for i in range(len(obj_list)):
        f.write(f"{i:5d} {obj_list[i]:.10f}\n")

print("end output", flush=True)
