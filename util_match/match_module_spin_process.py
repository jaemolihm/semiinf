from numpy import sin, cos
import numpy as np
from match_module_read_data import read_unkg, overlap_uk_to_ur, overlap_uk_to_ur_allrvec
dist_thr = 0.4
#################################
# some utility functions
def read_spnr(tbdict):
    # read seedname_spnr.bin file and leave only R=0 part
    path_seedname = tbdict['path'] + tbdict['seedname']
    nw = tbdict['nw']
    with open(path_seedname+'_spnr.bin', 'rb') as f:
        farr = np.fromfile(f, dtype=complex)
        spnr_temp = farr.reshape((nw, nw, tbdict['nrpts_all'], 3), order='F')
    tbdict['spnr_origin'] = spnr_temp[:,:,tbdict['ir_origin_all'],:].copy()
    tbdict['spnr_origin_backup'] = tbdict['spnr_origin'].copy()
    spnr_temp = None
    return

def cart2sph(xyz):
    xy = xyz[0]**2 + xyz[1]**2
    r = np.sqrt(xy + xyz[2]**2)
    theta = np.arctan2(np.sqrt(xy), xyz[2])
    phi = np.arctan2(xyz[1], xyz[0])
    return r, theta, phi

def spin_rotate(theta, phi, tbdict, iw0, iw1):
    # Return the unitary matrix to rotate given basis to z basis
    # calculate spin quantization axis in cartesian vector
    # We follow the convention of pw2wannier90.f90
    # n = (a,b,c)
    a = sin(theta) * cos(phi)
    b = sin(theta) * sin(phi)
    c = cos(theta)
    if 1 - abs(c) < 1E-2:
        # case |1-c| small(c close to 1): spin along z axis, assume no relative phase
        # This assumption is justified only if the original spinor phase of pw2wannier90.f90
        # is not changed.
        # if c negative, swap
        exp_phi = tbdict['spnr_origin'][iw0,iw1,0]
        if exp_phi.real > 0:
            phi = np.angle(exp_phi)
            fac0 = np.exp(1j*phi/2)
            fac1 = np.exp(-1j*phi/2)
        else:
            phi = np.angle(exp_phi) - np.pi
            fac0 = -np.exp(1j*phi/2)
            fac1 = np.exp(-1j*phi/2)
        if c > 0:
            return np.diag([fac1, fac0])
        else:
            return np.array([[0.0, fac0],[fac1, 0.0]], dtype=complex)
    else:
        # case |1-c| not small: spin not along z axis
        # calculate nz=(n cross z) = 1/sqrt(1-c*c)(b, -a, 0) and calculate coherence of s_nz
        # phase = 1j * np.exp(1j*delta)
        # where |0>=|up,n>, |1>=exp(1j*delta)*|dn,n>
        # |up,n> = 1/sqrt(2(1+c)) * (1+c, a+1j*b).T
        # |dn,n> = 1/sqrt(2(1-c)) * (-1+c, a+1j*b).T 
        # print(f"spin not along z axis for iw=({iw0:d},{iw1:d}): sx={a:.3f}, sy={b:.3f}, sz={c:.3f}")
        exp_phi = ( b * tbdict['spnr_origin'][iw0,iw1,0] - a * tbdict['spnr_origin'][iw0,iw1,1]) / np.sqrt(1-c*c) / 1j
        if exp_phi.real > 0:
            phi = np.angle(exp_phi)
            fac0 = np.exp(1j*phi/2)
            fac1 = np.exp(-1j*phi/2)
        else:
            phi = np.angle(exp_phi) - np.pi
            fac0 = -np.exp(1j*phi/2)
            fac1 = np.exp(-1j*phi/2)
        # u|up,z> = |0>, u|dn,z> = |1>
        uH = np.array([[   (1+c)/np.sqrt(2*(1+c)),   (-1+c)/np.sqrt(2*(1-c))],
                       [(a+1j*b)/np.sqrt(2*(1+c)), (a+1j*b)/np.sqrt(2*(1-c))]])
        uH[:,0] *= fac1
        uH[:,1] *= fac0
        u = uH.T.conjugate()
        return u
#################################

def tb_correct_spinor(tbbulk, tbslab, correct_sign):
    iw_head = tbslab['iw_head'] - tbslab['nw_match_add'] * 2
    iw_tail = tbslab['iw_tail'] + tbslab['nw_match_add'] * 2
    assert iw_head % 2 == 0
    assert iw_tail % 2 == 0
    read_unkg(tbbulk)
    read_unkg(tbslab)

    # if not spinor, swap and correct sign only.
    if (not tbbulk['isspinor']) or (not tbslab['isspinor']):
        tbslab['hk_spn'] = [x.copy() for x in tbslab['hk']]
        if correct_sign:
            set_manual_bulk_slab_pair(tbbulk, tbslab, iw_head, iw_tail)
            tb_correct_sign(tbbulk, tbslab)
        
        tbslab['hr_spn'] = overlap_uk_to_ur(tbslab['hk_spn'], tbslab, minussign=True)
        tbslab['hr_spn_all'] = overlap_uk_to_ur_allrvec(tbslab['hk_spn'], tbslab, minussign=True)
        return

    read_spnr(tbbulk)
    read_spnr(tbslab)
    spnr_b = tbbulk['spnr_origin']
    spnr_s = tbslab['spnr_origin']

    set_manual_bulk_slab_pair(tbbulk, tbslab, iw_head, iw_tail)

    # find spin pairs in bulk using spin operator value
    ind_pair_found = []
    pair_b = []
    for iw in range(0, tbbulk['nw']):
        if iw in ind_pair_found: continue
        spin_value = abs(spnr_b[:,iw,0])**2 + abs(spnr_b[:,iw,1])**2 + abs(spnr_b[:,iw,2])**2
        spin_value[iw] = 0
        iw_pair = np.argmax(spin_value)
        if iw_pair in ind_pair_found:
            print(f"tb_correct_spinor: problem in pair searching for iw={iw}, iw_pair={iw_pair}")
            raise ValueError
        pair_b += [(iw, iw_pair)]
        ind_pair_found += [iw, iw_pair]
        pair_dict = dict(pair_b + [(x[1],x[0]) for x in pair_b])

    ind_pair_found = []
    pair_s = []
    for iws in range(tbslab['nw']):
        if iws in ind_pair_found: continue
        # Find jws such that (iws - tbslab['iw_head']) % tbbulk['nw'] and
        # (jws - tbslab['iw_head']) % tbbulk['nw'] are pair_b pairs.
        # choose jws closest to iws
        iwb = (iws - tbslab['iw_head']) % tbbulk['nw']
        jwb_target = pair_dict[iwb]
        jws_cand = []
        for jws in range(tbslab['nw']):
            if (jws - tbslab['iw_head']) % tbbulk['nw'] == jwb_target:
                jws_cand.append(jws)
        jws_cand = np.array(jws_cand)
        jws_found = jws_cand[np.argmin(abs(jws_cand-iws))]
        pair_s += [(iws, jws_found)]
        ind_pair_found += [iws, jws_found]

    # get u_rotate from spin expectation value
    tbbulk['pair'] = pair_b
    tbslab['pair'] = pair_s

    # check if the WF center distance of spin pairs are close enough
    for p in pair_b:
        dist = np.linalg.norm(tbbulk['wfcenter'][:,p[0]] - tbbulk['wfcenter'][:,p[1]])
        if dist > dist_thr:
            print(f"WF center distance btw spin pairs is too far for bulk ({p[0]:d}, {p[1]:d}), distance={dist:.5f}")
            raise ValueError
    for p in pair_s:
        if (max(p) < iw_head) or (min(p) > iw_tail): continue
        dist = np.linalg.norm(tbslab['wfcenter'][:,p[0]] - tbslab['wfcenter'][:,p[1]])
        if dist > dist_thr:
            print(f"WF center distance btw spin pairs is too far for slab ({p[0]:d}, {p[1]:d}), distance={dist:.5f}")
            raise ValueError

    # if iw and iw+1 are not spin pair, print complaints
    for p in pair_b:
        spin_pairing = np.linalg.norm(spnr_b[p[1],p[0],:])**2
        if spin_pairing < 1.0:
            print(f"spin pairing is wrong for bulk ({p[0]:d}, {p[1]:d}), spin_pairing={spin_pairing:.3f}")
            raise ValueError
    for p in pair_s:
        if (max(p) < iw_head) or (min(p) > iw_tail): continue
        spin_pairing = np.linalg.norm(spnr_s[p[1],p[0],:])**2
        if spin_pairing < 1.0:
            print(f"spin pairing is wrong for slab ({p[0]:d}, {p[1]:d}), spin_pairing={spin_pairing:.3f}")
            raise ValueError

    u_rotate = np.eye(tbslab['nw'], dtype=complex)
    for p in pair_s:
        if (max(p) < iw_head) or (min(p) > iw_tail): continue
        iws0 = p[0]
        iws1 = p[1]
        iwb0 = (iws0 - tbslab['iw_head']) % tbbulk['nw']
        iwb1 = (iws1 - tbslab['iw_head']) % tbbulk['nw']

        spinvec_b = (spnr_b[iwb0,iwb0,:] - spnr_b[iwb1,iwb1,:]).real
        spinvec_s = (spnr_s[iws0,iws0,:] - spnr_s[iws1,iws1,:]).real
        r_b, theta_b, phi_b = cart2sph(spinvec_b)
        r_s, theta_s, phi_s = cart2sph(spinvec_s)

        U = ( spin_rotate(theta_b, phi_b, tbbulk, iwb0, iwb1)
            @ spin_rotate(theta_s, phi_s, tbslab, iws0, iws1).T.conjugate())
        u_rotate[iws0,iws0] = U[0,0]
        u_rotate[iws0,iws1] = U[0,1]
        u_rotate[iws1,iws0] = U[1,0]
        u_rotate[iws1,iws1] = U[1,1]

    # rotate operators using u_rotate
    u_rotate = np.asmatrix(u_rotate)
    for i in range(tbslab['sig'].shape[2]):
        tbslab['sig'][:,:,i] = np.dot(u_rotate, tbslab['sig'][:,:,i])
    tbslab['hk_spn'] = [u_rotate * x * u_rotate.H for x in tbslab['hk_spn']]
    for i in range(3):
        tbslab['spnr_origin'][:,:,i] = np.asarray(u_rotate * 
            np.asmatrix(tbslab['spnr_origin'][:,:,i]) * u_rotate.H)

    if correct_sign:
        tb_correct_phase(tbbulk, tbslab)

    tbslab['hr_spn'] = overlap_uk_to_ur(tbslab['hk_spn'], tbslab, minussign=True)
    tbslab['hr_spn_all'] = overlap_uk_to_ur_allrvec(tbslab['hk_spn'], tbslab, minussign=True)
    return

def tb_correct_sign(tbbulk, tbslab):
    isig_candidate = [0,1,2,3,4]
    iw_head = tbslab['iw_head'] - tbslab['nw_match_add'] * 2
    iw_tail = tbslab['iw_tail'] + tbslab['nw_match_add'] * 2

    sig_b = tbbulk['sig']
    sig_s = tbslab['sig']

    # use signature id 0,1,2,3,4 (average of WF in real space, cos(x), sin(x), cos(y), sin(y)) to get overlap matrix
    # use signature type with maximum absolute value
    # sin(z) cannot be used because Lz is different
    u_rotate = np.eye(tbslab['nw'], dtype=complex)
    
    for iws in range(iw_head, iw_tail):
        iwb = (iws - tbslab['iw_head']) % tbbulk['nw']

        sig_val = sig_s[iws, isig_candidate, :]
        isig_use, ispn_use = np.unravel_index(sig_val.argmax(), sig_val.shape)
        sig_iws = sig_s[iws, isig_use, ispn_use]
        sig_iwb = sig_b[iwb, isig_use, ispn_use]
        # calculate phase
        phase = 1 if (np.conjugate(sig_iws)*sig_iwb).real > 0 else -1
        u_rotate[iws, iws] = phase
        tbslab['sig'][iws,:,:] *= phase

    u_rotate = np.asmatrix(u_rotate)
    tbslab['hk_spn'] = [u_rotate.H * x * u_rotate for x in tbslab['hk_spn']]
    if tbbulk['isspinor']:
        for i in range(3):
            tbslab['spnr_origin'][:,:,i] = np.asarray(u_rotate.H * 
                np.asmatrix(tbslab['spnr_origin'][:,:,i]) * u_rotate)
    return


def tb_correct_phase(tbbulk, tbslab):
    isig_candidate = [0,1,2,3,4]
    iw_head = tbslab['iw_head'] - tbslab['nw_match_add'] * 2
    iw_tail = tbslab['iw_tail'] + tbslab['nw_match_add'] * 2

    sig_b = tbbulk['sig']
    sig_s = tbslab['sig']

    # use signature id 0,1,2,3,4 (average of WF in real space, cos(x), sin(x), cos(y), sin(y)) to get overlap matrix
    # use signature type with maximum absolute value
    # sin(z) cannot be used because Lz is different
    u_rotate = np.eye(tbslab['nw'], dtype=complex)
    
    for p in tbslab['pair']:
        if (max(p) < iw_head) or (min(p) > iw_tail): continue
        iws0 = p[0]
        iws1 = p[1]
        iwb0 = (iws0 - tbslab['iw_head']) % tbbulk['nw']
        iwb1 = (iws1 - tbslab['iw_head']) % tbbulk['nw']

        sig_iws = np.concatenate((sig_s[iws0, isig_candidate, :].flatten(), sig_s[iws1, isig_candidate, :].flatten()))
        sig_iwb = np.concatenate((sig_b[iwb0, isig_candidate, :].flatten(), sig_b[iwb1, isig_candidate, :].flatten()))
        # calculate phase
        phase = np.sum(sig_iwb * sig_iws.conjugate())
        if abs(phase) < 1E-4: print(f"abs(phase) used in phase matching might be too small: "
                                    f"({iws0:d},{iws1:d}), abs(phase) = {abs(phase):.3E}")
        phase = phase / abs(phase)
        u_rotate[iws0, iws0] = phase
        u_rotate[iws1, iws1] = phase
        tbslab['sig'][iws0,:,:] *= phase
        tbslab['sig'][iws1,:,:] *= phase

    u_rotate = np.asmatrix(u_rotate)
    tbslab['hk_spn'] = [u_rotate.H * x * u_rotate for x in tbslab['hk_spn']]
    if tbbulk['isspinor']:
        for i in range(3):
            tbslab['spnr_origin'][:,:,i] = np.asarray(u_rotate.H * 
                np.asmatrix(tbslab['spnr_origin'][:,:,i]) * u_rotate)
    return


def set_manual_bulk_slab_pair(tbbulk, tbslab, iw_head, iw_tail):
    path = tbbulk['path']
    tbslab['iws_to_iwb'] = {}
    for i in range(tbslab['nw']):
        tbslab['iws_to_iwb'][i] = (i - tbslab['iw_head']) % tbbulk['nw']
    if 'GeTe.5layer_1vac.conv9' in path:
        tbslab['iws_to_iwb'][112] = 10
        tbslab['iws_to_iwb'][113] = 11
        tbslab['iws_to_iwb'][114] = 8
        tbslab['iws_to_iwb'][115] = 9
        
        tbslab['iws_to_iwb'][121] = 21
        tbslab['iws_to_iwb'][125] = 17
        
        tbslab['iws_to_iwb'][73] = 21
        tbslab['iws_to_iwb'][77] = 17

        tbslab['iws_to_iwb'][160] = 10
        tbslab['iws_to_iwb'][161] = 11
        tbslab['iws_to_iwb'][162] = 8
        tbslab['iws_to_iwb'][163] = 9

        tbslab['iws_to_iwb'][169] = 21
        tbslab['iws_to_iwb'][173] = 17

        tbslab['iws_to_iwb'][176] = 26
        tbslab['iws_to_iwb'][177] = 27
        tbslab['iws_to_iwb'][178] = 24
        tbslab['iws_to_iwb'][179] = 25

        tbslab['iws_to_iwb'][217] = 21
        tbslab['iws_to_iwb'][221] = 17
    # elif path == '../main_data/GeTe/5layer_1vac.slab/conv_9_wrongspinx/':
    #     # tbslab['iws_to_iwb'][24] = 20
    #     # tbslab['iws_to_iwb'][28] = 16

    #     # tbslab['iws_to_iwb'][32] = 28
    #     # tbslab['iws_to_iwb'][33] = 29
    #     # tbslab['iws_to_iwb'][34] = 26
    #     # tbslab['iws_to_iwb'][35] = 27

    #     tbslab['iws_to_iwb'][50] = 44
    #     tbslab['iws_to_iwb'][52] = 42

    #     tbslab['iws_to_iwb'][72] = 20
    #     tbslab['iws_to_iwb'][76] = 16

    #     tbslab['iws_to_iwb'][98] = 44
    #     tbslab['iws_to_iwb'][100] = 42

    #     tbslab['iws_to_iwb'][112] = 10
    #     tbslab['iws_to_iwb'][113] = 11
    #     tbslab['iws_to_iwb'][114] = 8
    #     tbslab['iws_to_iwb'][115] = 9

    #     tbslab['iws_to_iwb'][120] = 20
    #     tbslab['iws_to_iwb'][124] = 16

    #     tbslab['iws_to_iwb'][146] = 44
    #     tbslab['iws_to_iwb'][148] = 42

    #     tbslab['iws_to_iwb'][160] = 10
    #     tbslab['iws_to_iwb'][161] = 11
    #     tbslab['iws_to_iwb'][162] = 8
    #     tbslab['iws_to_iwb'][163] = 9

    #     tbslab['iws_to_iwb'][168] = 20
    #     tbslab['iws_to_iwb'][172] = 16

    #     tbslab['iws_to_iwb'][176] = 26
    #     tbslab['iws_to_iwb'][177] = 27
    #     tbslab['iws_to_iwb'][178] = 24
    #     tbslab['iws_to_iwb'][179] = 25

    #     tbslab['iws_to_iwb'][192] = 43
    #     tbslab['iws_to_iwb'][193] = 44
    #     tbslab['iws_to_iwb'][194] = 40
    #     tbslab['iws_to_iwb'][195] = 41
    #     tbslab['iws_to_iwb'][196] = 42

    #     tbslab['iws_to_iwb'][208] = 10
    #     tbslab['iws_to_iwb'][209] = 11
    #     tbslab['iws_to_iwb'][210] = 8
    #     tbslab['iws_to_iwb'][211] = 9

    #     # tbslab['iws_to_iwb'][216] = 20
    #     # tbslab['iws_to_iwb'][220] = 16

    #     # tbslab['iws_to_iwb'][224] = 26
    #     # tbslab['iws_to_iwb'][225] = 27
    #     # tbslab['iws_to_iwb'][226] = 24
    #     # tbslab['iws_to_iwb'][227] = 25

    # elif (path in [f'../main_data/Bi2Se3/8layer_1vac.slab/conv_{conv}/'
    #                 for conv in range(7,10)]):
    #     tbslab['iws_to_iwb'][150] = 62
    #     tbslab['iws_to_iwb'][151] = 63
    #     tbslab['iws_to_iwb'][152] = 60
    #     tbslab['iws_to_iwb'][153] = 61

    #     tbslab['iws_to_iwb'][174] = 86
    #     tbslab['iws_to_iwb'][176] = 84

    #     tbslab['iws_to_iwb'][181] = 3
    #     tbslab['iws_to_iwb'][183] = 1
        
    #     tbslab['iws_to_iwb'][84] = 86
    #     tbslab['iws_to_iwb'][85] = 87
    #     tbslab['iws_to_iwb'][86] = 84
    #     tbslab['iws_to_iwb'][87] = 85
        
    # #     tbslab['iws_to_iwb'][73] = 21
    # #     tbslab['iws_to_iwb'][77] = 17

    # #     tbslab['iws_to_iwb'][160] = 10
    # #     tbslab['iws_to_iwb'][161] = 11
    # #     tbslab['iws_to_iwb'][162] = 8
    # #     tbslab['iws_to_iwb'][163] = 9

    # #     tbslab['iws_to_iwb'][169] = 21
    # #     tbslab['iws_to_iwb'][173] = 17

    # #     tbslab['iws_to_iwb'][176] = 26
    # #     tbslab['iws_to_iwb'][177] = 27
    # #     tbslab['iws_to_iwb'][178] = 24
    # #     tbslab['iws_to_iwb'][179] = 25


    final_indices = [-1] * tbslab['nw']
    for i in range(tbslab['nw']):
        # among j such that tbslab['iws_to_iwb'][j] = (i - tbslab['iw_head']) % tbbulk['nw'],
        # choose j closest to i: min(abs(j-i))
        j_cand = []
        for j in range(tbslab['nw']):
            if tbslab['iws_to_iwb'][j] == (i - tbslab['iw_head']) % tbbulk['nw']:
                j_cand.append(j)
        j_cand = np.array(j_cand)
        final_indices[i] = j_cand[np.argmin(abs(j_cand-i))]
    tbslab['final_indices'] = final_indices

    u_rotate = np.eye(tbslab['nw'], dtype=complex)[:,final_indices]
    tbslab['wfcenter'] = tbslab['wfcenter'][:,final_indices]
    tbslab['sig'] = tbslab['sig'][final_indices,:,:]
    tbslab['final_indices'] = final_indices

    u_rotate = np.asmatrix(u_rotate)
    # check unitarity of u_rotate
    if np.linalg.norm(u_rotate.H*u_rotate - np.eye(u_rotate.shape[0])) > 1E-10:
        print("tb_sort: u_rotate nor unitary")
        raise ValueError

    tbslab['hk_spn'] = [u_rotate.H * x * u_rotate for x in tbslab['hk']]
    tbslab['hr_spn'] = overlap_uk_to_ur(tbslab['hk_spn'], tbslab, minussign=True)
    if 'spnr_origin' in tbslab.keys():
        for i in range(3):
            tbslab['spnr_origin'][:,:,i] = np.asarray(u_rotate.H * 
                np.asmatrix(tbslab['spnr_origin'][:,:,i]) * u_rotate)
    return

# def tb_sort(tbbulk, tbslab, tboverlap, sort_indices=[]):
#     # # one might need to add some pairs manually in the future...
#     # path = tbbulk['path']
#     # if ('../main_data/GeTe/5layer_1vac.slab/conv_' in path):
#         # pairs_manual = []

#     # sanity check for sort_indices
#     if sorted([s[0] for s in sort_indices]) != sorted([s[1] for s in sort_indices]):
#         print("initial and final indices are not same as set")
#         raise ValueError

#     # do sort: for s in sort_indices: send WF_s[0] to WF_s[1]
#     sort_dict = dict(sort_indices)
#     final_indices = [-1] * tbslab['nw']
#     for iw in range(tbslab['nw']):
#         if iw in sort_dict.keys():
#             final_indices[iw] = sort_dict[iw]
#         else:
#             final_indices[iw] = iw
#
#     u_rotate = np.eye(tbslab['nw'], dtype=complex)[final_indices,:]
#     tbslab['wfcenter'] = tbslab['wfcenter'][:,final_indices]

#     u_rotate = np.asmatrix(u_rotate)
    
#     # check unitarity of u_rotate
#     if np.linalg.norm(u_rotate.H*u_rotate - np.eye(u_rotate.shape[0])) > 1E-10:
#         print("tb_sort: u_rotate nor unitary")
#         raise ValueError

#     tboverlap['uk_spn'] = [u_rotate * x for x in tboverlap['uk']]
#     tbslab['hk_spn'] = [u_rotate * x * u_rotate.H for x in tbslab['hk']]
#     tboverlap['ur_spn'] = overlap_uk_to_ur(tboverlap['uk_spn'], tbslab)
#     tbslab['hr_spn'] = overlap_uk_to_ur(tbslab['hk_spn'], tbslab, minussign=True)
#     if 'spnr_origin' in tbslab.keys():
#         for i in range(3):
#             tbslab['spnr_origin'][:,:,i] = np.asarray(u_rotate * 
#                 np.asmatrix(tbslab['spnr_origin'][:,:,i]) * u_rotate.H)
#     return

