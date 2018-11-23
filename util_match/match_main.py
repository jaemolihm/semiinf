import sys
from time import localtime, strftime
import numpy as np
import numpy.matlib
import scipy.linalg

from match_input import setup_input
from match_module_read_data import (tb_read_data, tb_sort_wcenter, tb_shift_fermi,
    tb_correct_data)
from match_module_spin_process import tb_correct_spinor
from match_module_ham_process import run_opt_diagonly
from match_module_ham_mismatch import (ham_mismatch_intra, ham_mismatch_inter,
                                       ham_maximal_hopping)

def write_hr_dat(tbslab, hr_type, postfix="hr_match"):
    hr_filehame = tbslab['path'] + tbslab['seedname'] + '_'+postfix+'.dat'
    assert len(tbslab[hr_type]) == tbslab['nrpts_all']
    assert tbslab[hr_type][0].shape == (tbslab['nw'], tbslab['nw'])

    with open(hr_filehame, 'w') as f:
        header = 'written by match_main.py at ' + strftime("%a, %d %b %Y %H:%M:%S +0000", localtime()) + '\n'

        f.write(header) # Date and time
        f.write(f"{tbslab['nw']:10d}\n")
        f.write(f"{tbslab['nrpts_all']:10d}")
        for ir in range(tbslab['nrpts_all']):
            if ir % 15 == 0: f.write("\n")
            f.write(f"{tbslab['ndegen_all'][ir]:5d}")
        f.write("\n")
        for ir in range(tbslab['nrpts_all']):
            for iw in range(tbslab['nw']):
                for jw in range(tbslab['nw']):
                    f.write(f"{tbslab['rvec_all'][0,ir]:5d}"
                            f"{tbslab['rvec_all'][1,ir]:5d}"
                            f"{tbslab['rvec_all'][2,ir]:5d}"
                            f"{jw+1:5d}{iw+1:5d}"
                            f"{tbslab[hr_type][ir][jw,iw].real:12.6f}"
                            f"{tbslab[hr_type][ir][jw,iw].imag:12.6f}")
                    f.write("\n")

################
# input parameters #
material, path, input_params = setup_input()
################################
# read binary files and initialize them
################################
print('begin initialization', flush=True)
tbbulk = tb_read_data(material=material, path=path, isbulk=True, input_params=input_params)
tbslab = tb_read_data(material=material, path=path, isbulk=False, input_params=input_params)
tb_shift_fermi(tbbulk, tbslab)
print('end initialization', flush=True)
################################
# for profiling
# import cProfile
# cProfile.run('run_opt_diagonly(tbbulk, tbslab, max_iter=200)')

# iden
ham_intra_max_iden, ham_intra_sqsum_iden = ham_mismatch_intra(tbslab['hr'], tbbulk, tbslab)
ham_inter_max_iden, ham_inter_sqsum_iden = ham_mismatch_inter(tbslab['hr'], tbbulk, tbslab)

# spn correction (minimal correction)
tb_correct_spinor(tbbulk, tbslab, correct_sign=True)
ham_intra_max_spn, ham_intra_sqsum_spn = ham_mismatch_intra(tbslab['hr_spn'], tbbulk, tbslab)
ham_inter_max_spn, ham_inter_sqsum_spn = ham_mismatch_inter(tbslab['hr_spn'], tbbulk, tbslab)

# hamiltonian correction
obj_list, is_converged = run_opt_diagonly(tbbulk, tbslab, max_iter=1000, init_spin_corrected=True)
ham_intra_max_ham, ham_intra_sqsum_ham = ham_mismatch_intra(tbslab['hr_ham'], tbbulk, tbslab)
ham_inter_max_ham, ham_inter_sqsum_ham = ham_mismatch_inter(tbslab['hr_ham'], tbbulk, tbslab)

ham_intra, ham_inter = ham_maximal_hopping(tbbulk, tbslab)

## printing
print("begin output", flush=True)

# write result summary
with open(path+'out_match.txt', 'w') as f:
    f.write("########################################################\n")
    if not is_converged: f.write("#####WARNING: convergence not reached during minimization#####\n\n")
    f.write("Calculation main_match.py (wannier_paper branch) at " 
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
    f.write(f"ham_intra_max_spn  = {ham_intra_max_spn:.10f}\n")
    f.write(f"ham_intra_max_ham  = {ham_intra_max_ham:.10f}\n")
    f.write("\n")
    f.write("===hamiltonian mismatch, intra-layer, sum of square (eV^2)===\n")
    f.write(f"ham_intra_sqsum_iden = {ham_intra_sqsum_iden:.10f}\n")
    f.write(f"ham_intra_sqsum_spn  = {ham_intra_sqsum_spn:.10f}\n")
    f.write(f"ham_intra_sqsum_ham  = {ham_intra_sqsum_ham:.10f}\n")
    f.write("\n")
    f.write("===hamiltonian mismatch, inter-layer, maximum error (eV)===\n")
    f.write(f"ham_inter_max_iden = {ham_inter_max_iden:.10f}\n")
    f.write(f"ham_inter_max_spn  = {ham_inter_max_spn:.10f}\n")
    f.write(f"ham_inter_max_ham  = {ham_inter_max_ham:.10f}\n")
    f.write("\n")
    f.write("===hamiltonian mismatch, inter-layer, sum of square (eV^2)===\n")
    f.write(f"ham_inter_sqsum_iden = {ham_inter_sqsum_iden:.10f}\n")
    f.write(f"ham_inter_sqsum_spn  = {ham_inter_sqsum_spn:.10f}\n")
    f.write(f"ham_inter_sqsum_ham  = {ham_inter_sqsum_ham:.10f}\n")
    f.write("########################################################\n")

# output minimization result in gnuplot format
with open(path+'out_match_minimization.txt', 'w') as f:
    f.write("# iter ham_mismatch(eV^2)\n")
    for i in range(len(obj_list)):
        f.write(f"{i:5d} {obj_list[i]:.10f}\n")

# save output hamiltonian as text file
write_hr_dat(tbslab, 'hr_all', postfix='hr_nocorr')
write_hr_dat(tbslab, 'hr_spn_all', postfix='hr_minimal')
write_hr_dat(tbslab, 'hr_ham_all', postfix='hr_ham')

print("end output", flush=True)
