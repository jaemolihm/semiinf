import os
import sys
import copy
import pickle
import numpy as np

def read_matlist(filename, dim1, dim2):
    with open(filename, 'rb') as f:
        farr = np.fromfile(f, dtype=complex)
        nlen = len(farr)//dim1//dim2
        arr_temp = farr.reshape((dim1, dim2, nlen), order='F')
        matlist = []
        for i in range(nlen):
            matlist.append(np.asmatrix(arr_temp[:,:,i]))
        return matlist

def read_array(filename, dims, dtype=np.float):
    with open(filename, 'rb') as f:
        array_out = np.fromfile(f, dtype=dtype).reshape(dims, order='F')
    return array_out

def tb_set_constants(tbdict, input_params):
    if tbdict['isbulk']:
        tbdict['seedname'] = input_params['seedname_b']
        tbdict['nbnd'] = input_params['nbnd_b']
        tbdict['nw'] = input_params['nw_b']
        try:
            tbdict['nw_exclude'] = input_params['nw_exclude_b']
        except (NameError, KeyError):
            pass
        tbdict['alat'] = input_params['alat_b']
        try:
          tbdict['atom_pos'] = input_params['atom_b']
        except (NameError, KeyError):
          pass
        tbdict['nklist'] = input_params['nklist_b']
        tbdict['iw_fermi'] = slice(0, input_params['nw_b'])
        tbdict['iw_match'] = slice(0, input_params['nw_b'])
        tbdict['norbital_list'] = input_params['norbital_list']
        assert sum(input_params['norbital_list']) == input_params['nw_b']
    else:
        tbdict['seedname'] = input_params['seedname_s']
        tbdict['nbnd'] = input_params['nbnd_s']
        tbdict['nw'] = input_params['nw_s']
        try:
            tbdict['nw_exclude'] = input_params['nw_exclude_s']
        except (NameError, KeyError):
            pass
        tbdict['alat'] = input_params['alat_s']
        try:
          tbdict['atom_pos'] = input_params['atom_s']
        except (NameError, KeyError):
          pass
        tbdict['nklist'] = input_params['nklist_s']
        tbdict['iw_head'] = input_params['iw_head']
        tbdict['iw_tail'] = input_params['iw_tail']
        tbdict['iw_fermi'] = slice(input_params['iw_head'], input_params['iw_tail'])
        tbdict['iw_match'] = slice(input_params['iw_head'], input_params['iw_tail'])
        
        iatm_add = input_params['iatm_add']
        tbdict['iatm_add'] =iatm_add
        tbdict['nw_match_add'] = sum(input_params['norbital_list'][0:iatm_add])
        assert sum(input_params['norbital_list'][0:iatm_add]) == sum(input_params['norbital_list'][-iatm_add:])
    
    tbdict['isspinor'] = input_params['isspinor']
    tbdict['nk'] = tbdict['nklist'][0] * tbdict['nklist'][1] * tbdict['nklist'][2]

def tb_read_data(material, path, isbulk, input_params):
    # print(f"begin reading {material}_{'bulk' if isbulk else 'slab'}"
    #       f" from {path}")
    
    # set constants
    tbdict = {}
    tbdict['material'] = material
    tbdict['path'] = path
    tbdict['isbulk'] = isbulk
    tb_set_constants(tbdict, input_params)

    # read binary files
    nk = tbdict['nk']
    nw = tbdict['nw']
    nbnd = tbdict['nbnd']
    path_seedname = tbdict['path'] + tbdict['seedname']

    # uloc_temp = read_matlist(path_seedname+'_u_matrix.bin', nw, nw)
    # udis_temp = read_matlist(path_seedname+'_u_matrix_opt.bin', nbnd, nw)

    # # if some bands are excluded from outer window in w90
    # if 'nw_exclude' in tbdict.keys():
    #     if not all([np.linalg.norm(x[-tbdict['nw_exclude']:,:])<1E-6 for x in udis_temp]):
    #         exit("Something's wrong with nw_exclude")
    #     for iks in range(len(udis_temp)):
    #         udis_temp[iks] = np.append(udis_temp[iks][-tbdict['nw_exclude']:,:],
    #                                 udis_temp[iks][0:-tbdict['nw_exclude'],:], axis=0)

    # tbdict['umat'] = [dismat * locmat for dismat, locmat 
    #                   in zip(udis_temp, uloc_temp)]
    # uloc_temp = None; udis_temp = None

    tbdict['hr'] = read_matlist(path_seedname+'_hr.bin', nw, nw)
    nrpts = len(tbdict['hr']) # this is changed after removing degenerate rvec points
    tbdict['rvec'] = read_array(path_seedname+'_irvec.bin', (3, nrpts), dtype=np.int32)
    tbdict['ndegen'] = read_array(path_seedname+'_ndegen.bin', (nrpts), dtype=np.int32)
    
    try:
        tbdict['wfcenter'] = read_array(path_seedname+'_wcenter.bin', (3, nw), dtype=np.float)
    except:
        tbdict['wfcenter'] = np.loadtxt(path_seedname+'_WF_center.txt', dtype=float)[:,1:4].T

    tbdict['kvec'] = read_array(path_seedname+'_kpt_latt.bin', (3, nk), dtype=np.float)   
    tbdict['ik_gamma'] = np.where(  (tbdict['kvec'][0,:]==0) 
                                  & (tbdict['kvec'][1,:]==0)
                                  & (tbdict['kvec'][2,:]==0))[0].item()

    # inverse fourier transform hr to hk
    try:
        f = open(path_seedname+'_hk.pkl', 'rb')
    except FileNotFoundError: # file does not exist
        tbdict['hk'] = []
        for ik in range(nk):
            tbdict['hk'].append(np.asmatrix(np.zeros((nw, nw), dtype=complex)))
            for ir in range(nrpts):
                rdotk = sum([r * k for r, k in zip(tbdict['rvec'][:,ir], tbdict['kvec'][:,ik])])
                tbdict['hk'][-1] += tbdict['hr'][ir] * np.exp(1j * 2 * np.pi * rdotk) / tbdict['ndegen'][ir]
        f = open(path_seedname+'_hk.pkl', 'wb')
        pickle.dump(tbdict['hk'], f)
        f.close()
    else: # file exists
        tbdict['hk'] = pickle.load(f)
        f.close()
    # check hermiticity
    for ik in range(nk):
        if np.linalg.norm(tbdict['hk'][ik] - tbdict['hk'][ik].H) > 1E-8:
            print(f"read_data: {material}_{'bulk' if isbulk else 'slab'}, " 
                  f"hk not hermitian at ik = {ik}")
            raise ValueError

    # remove degenerate rvec points
    tbdict['hr_all'] = copy.deepcopy(tbdict['hr'])
    tbdict['rvec_all'] = copy.deepcopy(tbdict['rvec'])
    tbdict['ndegen_all'] = copy.deepcopy(tbdict['ndegen'])
    tbdict['nrpts_all'] = len(tbdict['hr_all'])
    rvec_found = []
    irvec_dup = []
    for ir in range(nrpts):
        rvec_modulus = tuple(tbdict['rvec'][:,ir] % tbdict['nklist'])
        if rvec_modulus in rvec_found: # duplicate irvec, remove hr
            irvec_dup.append(ir)
        else:
            rvec_found.append(rvec_modulus)
    for ir in sorted(irvec_dup, reverse=True):
        del tbdict['hr'][ir];
        tbdict['rvec'] = np.delete(tbdict['rvec'], ir, 1)
    tbdict['ir_origin'] = np.where((tbdict['rvec'][0,:] == 0)
                                 & (tbdict['rvec'][1,:] == 0)
                                 & (tbdict['rvec'][2,:] == 0))[0].item()
    tbdict['nrpts'] = len(tbdict['hr'])
    tbdict['ndegen'] = None

    # print("end reading")
    return tbdict

# def overlap_read_data(tbbulk, tbslab, wholebulk=False):
#     # read binary files
#     tboverlapdict = {}

#     # read overlap matrix between Kohn-Sham eigenstates from file
#     path = tbbulk['path']
#     overlap_sb = read_matlist(path+'../overlap.bin', tbslab['nbnd'], tbbulk['nbnd'])
#     overlap = [np.conjugate(x.transpose()) for x in overlap_sb]

#     # find correspondence between ikb and iks k-points
#     ikb_to_iks = []
#     for ikb in range(tbbulk['nk']):
#         for iks in range(tbslab['nk']):
#             if np.linalg.norm(tbslab['kvec'][0:2,iks]
#                              - tbbulk['kvec'][0:2,ikb]) < 1E-6:
#                 ikb_to_iks.append(iks)
#                 continue
#     iks_to_ikb = {}
#     for iks in range(tbslab['nk']):
#         iks_to_ikb[iks] = tuple([ikb for ikb, iks_target in enumerate(ikb_to_iks)
#                                  if iks == iks_target])

#     # Apply wannierization unitary matrix to overlap matrix between 
#     # Kohn-Sham eigenstates to get overlap between Wannier functions
#     nw_b = tbbulk['nw']
#     if wholebulk:
#         rvecz = [-1, 0, 1]
#         uk_overlap = []
#         for iks in range(tbslab['nk']):
#             overlap_iks = np.zeros((len(rvecz)*nw_b, tbslab['nw']), dtype=complex)
#             for irz in range(len(rvecz)):
#                 for ikb in iks_to_ikb[iks]:
#                     phase = np.exp(1j*2*np.pi*tbbulk['kvec'][2,ikb] * rvecz[irz])
#                     overlap_iks[irz*nw_b:(irz+1)*nw_b] += (phase 
#                         * tbbulk['umat'][ikb].H * (overlap[ikb] * tbslab['umat'][iks]))
#             overlap_iks /= np.sqrt(tbbulk['nk']/tbslab['nk'])
#             # u, s, v = np.linalg.svd(overlap_iks, full_matrices=0)
#             # uk_overlap.append((np.asmatrix(u) * np.asmatrix(v)).H)
#             uk_overlap.append(np.asmatrix(overlap_iks).H)
#     else: # not wholebulk
#         uk_overlap = []
#         uk_overlap_h = []
#         uk_overlap_t = []
#         nw_add = tbslab['nw_match_add']
#         for iks in range(tbslab['nk']):
#             overlap_iks = np.zeros((tbbulk['nw'], tbslab['nw']), dtype=complex)
#             for ikb in iks_to_ikb[iks]:
#                 phase = np.exp(1j*2*np.pi*tbbulk['kvec'][2,ikb] * 0)
#                 overlap_iks += (phase * tbbulk['umat'][ikb].H * 
#                                         (overlap[ikb] * tbslab['umat'][iks]))
#             overlap_iks /= np.sqrt(tbbulk['nk']/tbslab['nk'])
#             uk_overlap.append(np.asmatrix(overlap_iks).H)
#             # upper slab
#             overlap_iks = np.zeros((nw_add, tbslab['nw']), dtype=complex)
#             for ikb in iks_to_ikb[iks]:
#                 phase = np.exp(1j*2*np.pi*tbbulk['kvec'][2,ikb] * 1)
#                 overlap_iks += (phase * (tbbulk['umat'][ikb][:,-nw_add:]).H * 
#                                         (overlap[ikb] * tbslab['umat'][iks]))
#             overlap_iks /= np.sqrt(tbbulk['nk']/tbslab['nk'])
#             uk_overlap_h.append(np.asmatrix(overlap_iks).H)
#             # lower slab
#             overlap_iks = np.zeros((nw_add, tbslab['nw']), dtype=complex)
#             for ikb in iks_to_ikb[iks]:
#                 phase = np.exp(1j*2*np.pi*tbbulk['kvec'][2,ikb] * -1)
#                 overlap_iks += (phase * (tbbulk['umat'][ikb][:,:nw_add]).H * 
#                                         (overlap[ikb] * tbslab['umat'][iks]))
#             overlap_iks /= np.sqrt(tbbulk['nk']/tbslab['nk'])
#             uk_overlap_t.append(np.asmatrix(overlap_iks).H)
        
#     # Fourier transform uk_overlap to ur_overlap
#     ur_overlap = []
#     for ir in range(tbslab['nrpts']):
#         ur_overlap.append(np.asmatrix(np.zeros(uk_overlap[0].shape, dtype=complex)))
#         for iks in range(tbslab['nk']):
#             rdotk = sum([r * k for r, k in zip(tbslab['rvec'][:,ir], tbslab['kvec'][:,iks])])
#             ur_overlap[-1] += uk_overlap[iks] * np.exp(1j * 2 * np.pi * rdotk)
#         ur_overlap[-1] /= tbslab['nk']

#     tboverlapdict['overlap'] = overlap
#     tboverlapdict['ikb_to_iks'] = ikb_to_iks
#     tboverlapdict['iks_to_ikb'] = iks_to_ikb
#     tboverlapdict['uk'] = uk_overlap
#     tboverlapdict['uk_head'] = uk_overlap_h
#     tboverlapdict['uk_tail'] = uk_overlap_t
#     tboverlapdict['ur'] = ur_overlap

    
#     if wholebulk:
#         raise ValueError # implementation not verified
#         tboverlapdict['iden'] = np.asmatrix(np.sign(tboverlapdict['ur'][tbslab['ir_origin']])
#                                             *(abs(tboverlapdict['ur'][tbslab['ir_origin']])>0.8  ))
#     else:
#         iden1 = np.zeros((tbslab['iw_head'], tbbulk['nw']), dtype=complex)
#         iden2 = np.eye(tbbulk['nw'], dtype=complex)
#         iden3 = np.zeros((tbslab['nw']-tbslab['iw_tail'], tbbulk['nw']), dtype=complex)
#         tboverlapdict['iden'] = np.asmatrix(np.concatenate((iden1, iden2, iden3), axis=0))

#     return tboverlapdict

def tb_shift_fermi(tb1, tb2):
    efermi_1 = np.average(np.diag(tb1['hr'][tb1['ir_origin']]
                                 [tb1['iw_fermi'], tb1['iw_fermi']]))
    efermi_2 = np.average(np.diag(tb2['hr'][tb2['ir_origin']]
                                 [tb2['iw_fermi'], tb2['iw_fermi']]))
    tb2['hr'][tb2['ir_origin']] += (efermi_1 - efermi_2) * np.eye(tb2['nw'])
    for hk_mat in tb2['hk']:
        hk_mat += (efermi_1 - efermi_2) * np.eye(tb2['nw'])

# def update_umat_s(tbbulk, tbslab, tboverlap, addpath):
#     tbslab_new = copy.deepcopy(tbslab)
#     tboverlap_new = copy.deepcopy(tboverlap)

#     # update tbslab
#     path_seedname = tbslab['path'] + addpath + tbslab['seedname']
#     uloc_temp = read_matlist(path_seedname+'_u_matrix.bin', tbslab['nw'], tbslab['nw'])
#     udis_temp = read_matlist(path_seedname+'_u_matrix_opt.bin', tbslab['nbnd'], tbslab['nw'])
    
#     # if some bands are excluded from outer window in w90
#     if 'nw_exclude' in tbdict.keys():
#         if not all([np.linalg.norm(x[-tbdict['nw_exclude']:,:])<1E-6 for x in udis_temp]):
#             exit("Something's wrong with nw_exclude")
#         for iks in range(len(udis_temp)):
#             udis_temp[iks] = np.append(udis_temp[iks][-tbdict['nw_exclude']:,:],
#                                     udis_temp[iks][0:-tbdict['nw_exclude'],:], axis=0)
    
#     tbslab_new['umat'] = [dismat * locmat for dismat, locmat 
#                       in zip(udis_temp, uloc_temp)]
#     uloc_temp = None; udis_temp = None


#     # update overlap
#     overlap = tboverlap['overlap']
#     uk_overlap = []
#     for iks in range(tbslab['nk']):
#         overlap_iks = np.zeros((tbbulk['nw'], tbslab['nw']), dtype=complex)
#         for ikb in tboverlap['iks_to_ikb'][iks]:
#             overlap_iks += tbbulk['umat'][ikb].H * (overlap[ikb] * tbslab_new['umat'][iks])
#         u, s, v = np.linalg.svd(overlap_iks, full_matrices=0)
#         uk_overlap.append((np.asmatrix(u) * np.asmatrix(v)).H)

#     ur_overlap = overlap_uk_to_ur(uk_overlap, tbslab)

#     tboverlap_new['uk'] = uk_overlap
#     tboverlap_new['ur'] = ur_overlap

#     return tbslab_new, tboverlap_new

def tb_sort_wcenter(tbbulk, tbslab):
    # sort wf indices using their center position to prevent swapping
    dist_thr = 0.2

    # try if iwb and iws=iwb+iw_head is a match
    wf_match_found_b = [False for iwb in range(tbbulk['nw'])]
    wf_match_found_s = [False for iwb in range(tbbulk['nw'])]
    ind_b_to_s = [iwb + tbslab['iw_head'] for iwb in range(tbbulk['nw'])]
    for iwb in range(tbbulk['nw']):
        dist = np.linalg.norm(tbbulk['wfcenter'][:,iwb]
                             - tbslab['wfcenter'][:,ind_b_to_s[iwb]])
        if dist < dist_thr:
            wf_match_found_b[iwb] = True
            wf_match_found_s[iwb] = True

    # if previous try fails, find other match
    ind_b_to_s_sorted = [x for x in ind_b_to_s]
    for iwb in range(tbbulk['nw']):
        if wf_match_found_b[iwb]: continue # match is already found
        for jwb in range(tbbulk['nw']):
            if wf_match_found_s[jwb]: continue # match is already found
            dist = np.linalg.norm(tbbulk['wfcenter'][:,iwb]
                                 - tbslab['wfcenter'][:,ind_b_to_s[jwb]])
            if dist < dist_thr:
                wf_match_found_b[iwb] = True
                wf_match_found_s[jwb] = True
                ind_b_to_s_sorted[iwb] = ind_b_to_s[jwb]
                break # goto next iwb
      
        if not wf_match_found_b[iwb]: # match is not found
            print(f"Error in tb_sort_wcenter: Match between bulk and slab "
                  f"is not found for iwb = {iwb}")
            raise ValueError

    # change umat, hr, wfcenter of bulk according to ind_b_to_s_sorted
    ind_s_new = [i for i in range(tbslab['nw'])]
    for iwb, iws in enumerate(ind_b_to_s_sorted):
        ind_s_new[ind_b_to_s[iwb]] = iws
    # tbslab['umat'] = [umat[:,ind_s_new] for umat in tbslab['umat']]
    tbslab['hr'] = [hrmat[ind_s_new,:][:,ind_s_new] for hrmat in tbslab['hr']]
    tbslab['hr_all'] = [hrmat[ind_s_new,:][:,ind_s_new] for hrmat in tbslab['hr_all']]
    tbslab['wfcenter'] = tbslab['wfcenter'][:,ind_s_new]

    #update hk according to updated hr
    # inverse fourier transform hr to hk
    tbslab['hk'] = []
    for ik in range(tbslab['nk']):
        tbslab['hk'].append(np.asmatrix(np.zeros((tbslab['nw'], tbslab['nw']), dtype=complex)))
        for ir in range(tbslab['nrpts_all']):
            rdotk = sum([r * k for r, k in zip(tbslab['rvec_all'][:,ir], tbslab['kvec'][:,ik])])
            tbslab['hk'][-1] += tbslab['hr_all'][ir] * np.exp(1j * 2 * np.pi * rdotk) / tbslab['ndegen_all'][ir]
    os.remove(tbslab['path'] + tbslab['seedname']+'_hk.pkl') # do not allow reuse
    # check hermiticity
    for ik in range(tbslab['nk']):
        norm = np.linalg.norm(tbslab['hk'][ik] - tbslab['hk'][ik].H)
        if norm > 1E-8:
            print(f"read_data: {tbslab['material']}_'slab', " 
                  f"hk not hermitian at ik = {ik}: {norm}")
            raise ValueError


def overlap_uk_to_ur(mat_k_list, tbslab, minussign=False):
    # Fourier transform mat_k_list to mat_r_list
    # minussign=True for hamiltonian
    # minussign=False for overlap
    mat_r_list = []
    for ir in range(tbslab['nrpts']):
        mat_r_list.append(np.asmatrix(np.zeros(mat_k_list[0].shape, dtype=complex)))
        for iks in range(tbslab['nk']):
            rdotk = sum([r * k for r, k in zip(tbslab['rvec'][:,ir], tbslab['kvec'][:,iks])])
            if minussign:
                mat_r_list[-1] += mat_k_list[iks] * np.exp(-1j * 2 * np.pi * rdotk)
            else:
                mat_r_list[-1] += mat_k_list[iks] * np.exp(1j * 2 * np.pi * rdotk)
        mat_r_list[-1] /= tbslab['nk']
    return mat_r_list

def overlap_uk_to_ur_allrvec(mat_k_list, tbslab, minussign=False):
    # Fourier transform mat_k_list to mat_r_list
    # minussign=True for hamiltonian
    # minussign=False for overlap
    mat_r_list = []
    for ir in range(tbslab['nrpts_all']):
        mat_r_list.append(np.asmatrix(np.zeros(mat_k_list[0].shape, dtype=complex)))
        for iks in range(tbslab['nk']):
            rdotk = sum([r * k for r, k in zip(tbslab['rvec_all'][:,ir], tbslab['kvec'][:,iks])])
            if minussign:
                mat_r_list[-1] += mat_k_list[iks] * np.exp(-1j * 2 * np.pi * rdotk)
            else:
                mat_r_list[-1] += mat_k_list[iks] * np.exp(1j * 2 * np.pi * rdotk)
        mat_r_list[-1] /= tbslab['nk']
    return mat_r_list

# def read_unkg(tbdict):
#     N_SIG = 7
#     # read and calculate signature
#     nbnd = tbdict['nbnd']
#     nw = tbdict['nw']
#     ikgamma = tbdict['ik_gamma']
#     isspinor = tbdict['isspinor']

#     g_abc = [[ 0, 0, 0],
#              [ 1, 0, 0],
#              [-1, 0, 0],
#              [ 0, 1, 0],
#              [ 0,-1, 0],
#              [ 0, 0, 1],
#              [ 0, 0,-1]]
#     g_abc_to_isig = {}
#     for isig in range(N_SIG):
#         g_abc_to_isig[tuple(g_abc[isig])] = isig

#     if isspinor: unkg = np.empty((nbnd, N_SIG, 2), dtype=complex)
#     if not isspinor: unkg = np.empty((nbnd, N_SIG, 1), dtype=complex)

#     is_spin_up = True # used only if isspinor==True
#     with open(tbdict['path']+'../'+tbdict['seedname']+'.unkg', 'r') as f:
#         lines = f.readlines()
#         for line in lines[1:]:
#             line2 = line.split()
#             g_abc_tuple = tuple([int(x) for x in line2[2:5]])
#             if g_abc_tuple in g_abc_to_isig.keys():
#                 isig = g_abc_to_isig[g_abc_tuple]
#                 if isspinor and is_spin_up:
#                     unkg[int(line2[0])-1, isig, 0] = float(line2[-2]) + 1j*float(line2[-1])
#                     is_spin_up = False
#                 elif isspinor:
#                     unkg[int(line2[0])-1, isig, 1] = float(line2[-2]) + 1j*float(line2[-1])
#                     is_spin_up = True
#                 else:
#                     unkg[int(line2[0])-1, isig, 0] = float(line2[-2]) + 1j*float(line2[-1])

#     alat_inv = np.linalg.inv(tbdict['alat'])
#     if isspinor: sig = np.zeros((nw, N_SIG, 2), dtype=complex)
#     if not isspinor: sig = np.zeros((nw, N_SIG, 1), dtype=complex)
#     for iw in range(nw):
#         for ib in range(nbnd):
#             sig[iw, :, :] += tbdict['umat'][tbdict['ik_gamma']][ib,iw] * unkg[ib, :, :]
    
#     tbdict['g_dot_wffrac'] = np.zeros((nw,N_SIG))
#     for iw in range(nw):
#         wf_frac = np.dot(alat_inv, tbdict['wfcenter'][:,iw])
#         for isig in range(N_SIG):
#             tbdict['g_dot_wffrac'][iw,isig] = sum(np.multiply(g_abc[isig], wf_frac))
#             phase_factor = np.exp(-1j*2*np.pi*sum(np.multiply(g_abc[isig], wf_frac)))
#             sig[iw,isig] *= phase_factor
    
#     if isspinor: sig_real = np.zeros((nw, N_SIG, 2), dtype=complex)
#     if not isspinor: sig_real = np.zeros((nw, N_SIG, 1), dtype=complex)
#     sig_real[:,0,:] = sig[:,0,:]
#     for isig in range(1, N_SIG//2+1):
#         # cos(x), cos(y), cos(z)
#         sig_real[:,2*isig-1, :] = (sig[:, 2*isig-1, :] + sig[:, 2*isig, :])/2
#         # sin(x), sin(y), sin(z)
#         sig_real[:,2*isig-1, :] = (sig[:, 2*isig-1, :] - sig[:, 2*isig, :])/(2*1j)

#     tbdict['unkg'] = unkg
#     tbdict['sig'] = sig_real


def tb_correct_data(tbbulk, tbslab):
    return
    # ad hoc correction for some materials: correct sign or WF swap
    path = tbbulk['path']
    do_parity = False
    u = None
    if ('../main_data/GeTe/5layer_1vac.relax/conv_' in path):
        do_parity = True
        u = np.asmatrix(np.eye(tbslab['nw'], dtype=complex))
        for i in [-22,10,11,17,26,27,42,43,90,91,69,75]:
            ii = i + tbslab['iw_head']
            u[ii,ii] = -1.0
        for i,j in [(-24,-22),(-23,-21),(56,58),(57,59),(72,74),(73,75),(88,90),(89,91)]:
            ii = i + tbslab['iw_head']
            jj = j + tbslab['iw_head']
            u[ii,ii] = 0.0
            u[jj,jj] = 0.0
            u[ii,jj] = -1.0
            u[jj,ii] = 1.0
            temp = tbslab['wfcenter'][:,ii].copy()
            tbslab['wfcenter'][:,ii] = tbslab['wfcenter'][:,jj].copy()
            tbslab['wfcenter'][:,jj] = temp.copy()
        for i,j in [(65,69),(-31,-27)]:
            ii = i + tbslab['iw_head']
            jj = j + tbslab['iw_head']
            u[ii,ii] = 0.0
            u[jj,jj] = 0.0
            u[ii,jj] = 1.0
            u[jj,ii] = -1.0
            temp = tbslab['wfcenter'][:,ii].copy()
            tbslab['wfcenter'][:,ii] = tbslab['wfcenter'][:,jj].copy()
            tbslab['wfcenter'][:,jj] = temp.copy()
    # if ('../main_data/GeTe/5layer_1vac.slab.old/conv_' in path):
    #     do_parity = True
    #     u = np.asmatrix(np.eye(tbslab['nw'], dtype=complex))
    #     for i in [10,11,17,56,57]:
    #         ii = i + tbslab['iw_head']
    #         u[ii,ii] = -1.0
    if ('../main_data/GeTe/5layer_1vac.slab/conv_' in path):
        do_parity = True
        u = np.asmatrix(np.eye(tbslab['nw'], dtype=complex))
        for i in [-27,10,11,17,56,57,69,72,73]:
            ii = i + tbslab['iw_head']
            u[ii,ii] = -1.0
        for i,j in [(-31,-27)]:
            ii = i + tbslab['iw_head']
            jj = j + tbslab['iw_head']
            u[ii,ii] = 0.0
            u[jj,jj] = 0.0
            u[ii,jj] = 1.0
            u[jj,ii] = -1.0
            temp = tbslab['wfcenter'][:,ii].copy()
            tbslab['wfcenter'][:,ii] = tbslab['wfcenter'][:,jj].copy()
            tbslab['wfcenter'][:,jj] = temp.copy()
    if ('../main_data/GeTe/7layer_1vac.relax/conv_' in path):
        do_parity = True
        u = np.asmatrix(np.eye(tbslab['nw'], dtype=complex))
        for i in [0,1,2,3,12,13,16,17,18,19,32,33]:
            ii = i + tbslab['iw_head']
            u[ii,ii] = -1.0
    if '../main_data/Bi2Se3/8layer_1vac.relax/conv_7' in path:
        do_parity = True
        u = np.asmatrix(np.eye(tbslab['nw'], dtype=complex))
        for i in [0,24,30,54,60,-90,-36,-30,-6,90,144]:
            ii = i + tbslab['iw_head']
            jj = ii + 2
            u[ii,ii] = 0.0
            u[jj,jj] = 0.0
            u[ii,jj] = -1.0
            u[jj,ii] = 1.0
            temp = tbslab['wfcenter'][:,ii].copy()
            tbslab['wfcenter'][:,ii] = tbslab['wfcenter'][:,jj].copy()
            tbslab['wfcenter'][:,jj] = temp.copy()
    if '../main_data/Bi2Se3/8layer_1vac.slab/conv_7' in path:
        do_parity = True
        u = np.asmatrix(np.eye(tbslab['nw'], dtype=complex))
        for i,j in [(84,86),(62,60),(63,61)]:
            ii = i + tbslab['iw_head']
            jj = j + tbslab['iw_head']
            u[ii,ii] = 0.0
            u[jj,jj] = 0.0
            u[ii,jj] = -1.0
            u[jj,ii] = 1.0
            temp = tbslab['wfcenter'][:,ii].copy()
            tbslab['wfcenter'][:,ii] = tbslab['wfcenter'][:,jj].copy()
            tbslab['wfcenter'][:,jj] = temp.copy()
        for i in [60, 62, 93, 122, 146, 147, -3, -4, -28, -33, -64, -87]:
            ii = i + tbslab['iw_head']
            u[:,ii] *= -1
    if do_parity:
        tbslab['umat'] = [x*u for x in tbslab['umat']]
        tbslab['hr'] = [u.H*x*u for x in tbslab['hr']]
        tbslab['hr_all'] = [u.H*x*u for x in tbslab['hr_all']]
        tbslab['hk'] = [u.H*x*u for x in tbslab['hk']]
    return
