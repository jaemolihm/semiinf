from time import localtime, strftime
import numpy as np
import numpy.matlib
import scipy.linalg
import pickle

def ham_mismatch_intra(hr_s, tbbulk, tbslab):
    diff_max_list = []
    diff_norm = 0.0
    for ir_s in range(tbslab['nrpts']):
        ir_b = np.where((tbbulk['rvec'][0,:] == tbslab['rvec'][0,ir_s])
                      & (tbbulk['rvec'][1,:] == tbslab['rvec'][1,ir_s])
                      & (tbbulk['rvec'][2,:] == 0))[0].item()
        hr_diff = hr_s[ir_s][tbslab['iw_match'],tbslab['iw_match']] - tbbulk['hr'][ir_b]
        diff_max_list.append(abs(hr_diff).max())
        diff_norm += np.linalg.norm(hr_diff) ** 2
    return max(diff_max_list), diff_norm

def ham_mismatch_inter(hr_s, tbbulk, tbslab):
    diff_max_list = []
    diff_norm = 0.0
    above = True # False: below
    # above = False: not sure this works correctly.
    if above:
        iw_intra = slice(tbslab['iw_head']-tbbulk['nw'], tbslab['iw_head'])
    else:
        iw_intra = slice(tbslab['iw_tail'], tbslab['iw_tail']+tbbulk['nw'])
    for ir_s in range(tbslab['nrpts']):
        ir_b = np.where((tbbulk['rvec'][0,:] == tbslab['rvec'][0,ir_s])
                      & (tbbulk['rvec'][1,:] == tbslab['rvec'][1,ir_s])
                      & (tbbulk['rvec'][2,:] == (1 if above else -1)))[0].item()
        hr_diff = hr_s[ir_s][tbslab['iw_match'],iw_intra] - tbbulk['hr'][ir_b]
        diff_max_list.append(abs(hr_diff).max())
        diff_norm += np.linalg.norm(hr_diff) ** 2
    return max(diff_max_list), diff_norm

def ham_maximal_hopping(tbbulk, tbslab):
    intra_max_list = []
    inter_max_list = []
    for ir_b in range(tbbulk['nrpts']):
        if tbbulk['rvec'][2,ir_b] == 0:
            hr_temp = tbbulk['hr'][ir_b]
            hr_temp = hr_temp - np.diag(np.diag(hr_temp))
            intra_max_list.append(abs(hr_temp).max())
        elif tbbulk['rvec'][2,ir_b] == 1: 
            inter_max_list.append(abs(tbbulk['hr'][ir_b]).max())
    return max(intra_max_list), max(inter_max_list)


def write_hr_dat(tbslab):
    hr_filehame = tbslab['path'] + tbslab['seedname'] + '_hr_match.dat'
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
                            f"{tbslab['hr_ham_all'][ir][jw,iw].real:12.6f}"
                            f"{tbslab['hr_ham_all'][ir][jw,iw].imag:12.6f}")
                    f.write("\n")


