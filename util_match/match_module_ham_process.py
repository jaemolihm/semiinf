import numpy as np
import numpy.matlib
import scipy.linalg
import pickle

from match_module_read_data import overlap_uk_to_ur, overlap_uk_to_ur_allrvec

# IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 50, NO. 3, MARCH 2002    
# "Optimization Algorithms Exploiting Unitary Constraints"
# Modified Steepest Descent on the Complex Grassmann Manifold 
# in subsection V-A, p.641
# modified to incorporate multiple co-optimization of multiple matrices
match_dist_threshold = 0.2 # angstrom
is_setup_H1_H2 = False
def phase_fix(U, rng1, rng2):
    # fix sum_ik U[ik][ind1,ind2] to be positive real number
    # by multiplying constant phase factor
    phase_fix_val = sum([sum(np.diag(umat[rng1,rng2])) for umat in U])
    phase_fix = abs(phase_fix_val) / phase_fix_val
    for umat in U: umat *= phase_fix

# def gradient(H1, H2, U, locality_args):
#     Gamma_locality = locality_penalty_grad(U, locality_args)
#     Z = []
#     for i in range(len(U)):
#         gamma = 2*(H2[i]*U[i])*(U[i].H*H2[i]*U[i] - H1[i]) \
#                 + locality_args['sigma'] * Gamma_locality[i]
#         gamma = np.asmatrix(gamma)
#         Z.append( U[i]*(gamma.H*U[i]) - gamma )
#     return Z

def objective(H1, H2, U):
    value = 0.0
    for i in range(len(H1)):
        uh2u = U[i].H * H2[i] * U[i]
        value += -2 * np.einsum('ij,ji->', H1[i], uh2u).real
        # the below contribution is added in obj_offset
        # this is possible only if U[i] is a square, unitary matrix
        # value += np.linalg.norm(uh2u)**2
    return value

# def project_to_unitary(U):
#     for i in range(len(U)):
#         u, s, v = np.linalg.svd(U[i], full_matrices=False)
#         U[i] = np.asmatrix(np.dot(u, v))

def block_diag_expm(mat, rng_list):
    exp_mat = np.zeros(mat.shape, dtype=complex)
    for rng in rng_list:
        exp_mat[rng,rng] = scipy.linalg.expm(mat[rng,rng])
    return exp_mat


def run_opt_diagonly(tbbulk, tbslab, do_unit_orbital=True, max_iter=None,
                     init_spin_corrected=True):
    nk_s = tbslab['nk']
    nw_s = tbslab['nw']
    nk_b = tbbulk['nk']
    nw_b = tbbulk['nw']
    verbose = True
    if max_iter is None: # set to default value
        max_iter = 1000
        verbose = True

    # setup index match, H1, H2, use set values if done before.
    set_index_match(tbbulk, tbslab)
    H1, H2 = setup_H1_and_H2(tbbulk, tbslab, init_spin_corrected=init_spin_corrected)

    U = [np.asmatrix(np.eye(tbslab['nw_match']), dtype=complex) for i in range(nk_s)]

    # block diagonality constraint on the shape of U
    # assume tbdict['nw_match_add'] = sum(norbital_list[0:4]) has not changed
    iatm_add = tbslab['iatm_add']
    if do_unit_orbital:
        nunit = (tbbulk['norbital_list']
            + tbbulk['norbital_list'][-iatm_add:] # above iw_head
            + tbbulk['norbital_list'][0:iatm_add]) # below iw_tail
    elif tbbulk['is_spinor']:
        nunit = [2] * (nw_b+tbslab['nw_match_add'])//2
    else:
        nunit = [1] * (nw_b+tbslab['nw_match_add'])
    nunit_cumsum = np.append(0, np.cumsum(np.array(nunit)))
    if nunit_cumsum[-1] != H1[0].shape[0]: raise ValueError
    rng_list = []
    for i in range(len(nunit_cumsum)-1):
        rng_list += [slice(nunit_cumsum[i], nunit_cumsum[i+1])]

        
    mu = 2**-10
    thr = 1E-5
    norm_G = thr+1
    niter = 0
    obj_offset = sum([np.linalg.norm(hmat)**2 for hmat in H1]).real
    obj_offset += sum([np.linalg.norm(hmat)**2 for hmat in H2]).real
    obj_list = [objective(H1, H2, U)]
    is_converged = True
    if verbose: print("begin iterative minimization", flush=True)
    is_converged = False
    for niter in range(max_iter):
        gamma = np.asmatrix(np.zeros(U[0].shape, dtype=complex))
        for i in range(len(H1)):
            gamma += -2*H2[i]*U[i]*H1[i]
        G = np.asmatrix(np.zeros(U[0].shape, dtype=complex))
        for rng in rng_list: G[rng,rng] = ( gamma[rng,rng] * U[0][rng,rng].H
                                          - U[0][rng,rng] * gamma[rng,rng].H )
        norm_G = 0.5*np.sum(G*G.H).real
        # step 5: update gradient step size mu if too small
        P = block_diag_expm(-mu*G, rng_list)
        Q = block_diag_expm(-2*mu*G, rng_list)
        temp_mat = Q*U[0]
        QU = [temp_mat for i in range(len(U))]
        dobj_q = obj_list[-1] - objective(H1, H2, QU)
        while dobj_q >= mu * norm_G:
            mu *= 2
    #        print("mu = ", mu)
            Q = block_diag_expm(-2*mu*G, rng_list)
            temp_mat = Q*U[0]
            QU = [temp_mat for i in range(len(U))]
            dobj_q = obj_list[-1] - objective(H1, H2, QU)
        #step 7
        temp_mat = P*U[0]
        PU = [temp_mat for i in range(len(U))]
        dobj_p = obj_list[-1] - objective(H1, H2, PU)
        while (dobj_p < mu/2 * norm_G):
            mu *= 0.5
    #        print("mu = ", mu)
            P = block_diag_expm(-mu*G, rng_list)
            temp_mat = P*U[0]
            PU = [temp_mat for i in range(len(U))]
            dobj_p = obj_list[-1] - objective(H1, H2, PU)
        # setp 7: set U, goto step 2
        U = PU
        obj_list.append(objective(H1, H2, U))
        if mu < 1E-30: break
        if norm_G < thr:
            print(f"niter={niter}, norm_G={norm_G:.1E}, "
                  f"dobj_q={dobj_q:.1E}, dobj_p={dobj_p:.1E}", flush=True)
            is_converged = True
            break
        if niter%100 == 0 and verbose:
            print(f"niter={niter}, norm_G={norm_G:.1E}, "
                  f"dobj_q={dobj_q:.1E}, dobj_p={dobj_p:.1E}", flush=True)
    if verbose: print("end iterative minimization", flush=True)
    obj_list = [x + obj_offset for x in obj_list]
  

    sl_c = tbslab['iw_match']
    sl_h = slice(tbslab['iw_head']-tbslab['nw_match_add'], tbslab['iw_head'])
    sl_t = slice(tbslab['iw_tail'], tbslab['iw_tail']+tbslab['nw_match_add'])
    sl_Uc = slice(0, nw_b)
    sl_Uh = slice(nw_b, nw_b + tbslab['nw_match_add'])
    sl_Ut = slice(nw_b + tbslab['nw_match_add'], nw_b + 2*tbslab['nw_match_add'])
    u_rotate = np.eye(tbslab['nw'], dtype=complex)
    u_rotate[sl_c, sl_c] = U[0][sl_Uc, sl_Uc]
    u_rotate[sl_h, sl_h] = U[0][sl_Uh, sl_Uh]
    u_rotate[sl_t, sl_t] = U[0][sl_Ut, sl_Ut]
    u_rotate = np.asmatrix(u_rotate)

    if init_spin_corrected:
        tbslab['hk_ham'] = [u_rotate.H * x * u_rotate for x in tbslab['hk_spn']]
    else:
        tbslab['hk_ham'] = [u_rotate.H * x * u_rotate for x in tbslab['hk']]
    tbslab['hr_ham'] = overlap_uk_to_ur(tbslab['hk_ham'], tbslab, minussign=True)
    tbslab['hr_ham_all'] = overlap_uk_to_ur_allrvec(tbslab['hk_ham'], tbslab, minussign=True)
    return obj_list, is_converged


def set_index_match(tbbulk, tbslab):
    nw_b = tbbulk['nw']
    nw_match = nw_b + 2*tbslab['nw_match_add']
    ind_slab_match = list(range(tbslab['iw_head'], tbslab['iw_tail']))
    ind_slab_match += list(range(tbslab['iw_head']-tbslab['nw_match_add'], tbslab['iw_head']))
    ind_slab_match += list(range(tbslab['iw_tail'], tbslab['iw_tail']+tbslab['nw_match_add']))

    match_found = [True] * nw_b + [False] * (nw_match-nw_b)

    ind_layer_match = [0] * nw_b + [-999] * (nw_match-nw_b)
    ind_bulk_match = list(range(0, nw_b)) + [-999] * (nw_match-nw_b)


    # Now, find bulk WFs corresponding to ind_slab_match[nw_b:]
    # As a first guess, if ind_slab_match(j) +- num_wann_bulk is in ind_b_to_s,
    # try ind_bulk_match(j)=i & ilayer = -1 or +1
    for imatch in range(nw_b, nw_match):
        for iwb in range(0, nw_b):
            if not (ind_slab_match[imatch] + nw_b == ind_slab_match[iwb]
                    or ind_slab_match[imatch] - nw_b == ind_slab_match[iwb]): continue
            for ilayer in (-1,1):
                dist = np.linalg.norm(tbbulk['wfcenter'][:,iwb]
                                    - tbslab['wfcenter'][:,ind_slab_match[imatch]]
                                    + tbbulk['alat'][:,2] * ilayer)
                if dist < match_dist_threshold:
                    # check for duplicate
                    is_duplicate = False
                    for jmatch in range(nw_match):
                        if (iwb==ind_bulk_match[jmatch] 
                            and ilayer==ind_layer_match[jmatch]):
                            is_duplicate = True
                            break # break duplicate test
                    if not is_duplicate: # if not duplicate, match is found
                        match_found[imatch] = True
                        ind_bulk_match[imatch] = iwb
                        ind_layer_match[imatch] = ilayer
                        break # exit search for ilayer
            if match_found[imatch]: break # exit search for iwb

    # If previous guess was not successful, try all bulk wfs, with ilayer = -2 ~ 2
    for imatch in range(nw_b, nw_match):
        if match_found[imatch]: continue # match found in previous guess
        for iwb in range(0, nw_b):
            for ilayer in range(-2,3):
                dist = np.linalg.norm(tbbulk['wfcenter'][:,iwb]
                                    - tbslab['wfcenter'][:,ind_slab_match[imatch]]
                                    + tbbulk['alat'][:,2] * ilayer)
                if dist < match_dist_threshold:
                    # check for duplicate
                    is_duplicate = False
                    for jmatch in range(nw_match):
                        if (iwb==ind_bulk_match[jmatch] 
                            and ilayer==ind_layer_match[jmatch]):
                            is_duplicate = True
                            break # break duplicate test
                    if not is_duplicate: # if not duplicate, match is found
                        match_found[imatch] = True
                        ind_bulk_match[imatch] = iwb
                        ind_layer_match[imatch] = ilayer
                        break # exit search for ilayer
            if match_found[imatch]: break # exit search for iwb

    if not all(match_found):
        print("Error in module_ham_process set_index_match: "
              "Match between bulk and slab.\n"
              "Try increasing match_dist_threshold.")
        raise ValueError

    tbslab['nw_match'] = nw_match
    tbslab['ind_slab_match'] = ind_slab_match
    tbbulk['ind_bulk_match'] = ind_bulk_match
    tbbulk['ind_layer_match'] = ind_layer_match

def setup_H1_and_H2(tbbulk, tbslab, init_spin_corrected):
    nw_match = tbslab['nw_match']
    ind_slab_match = tbslab['ind_slab_match']
    ind_bulk_match = tbbulk['ind_bulk_match']
    ind_layer_match = tbbulk['ind_layer_match']

    hr_fourier = {}
    for ir in range(tbbulk['nrpts_all']):
        rvec_z = tbbulk['rvec_all'][2,ir]
        if rvec_z not in hr_fourier.keys():
            hr_fourier[rvec_z] = [np.zeros((tbbulk['nw'],tbbulk['nw']), dtype=complex) for ik in range(tbslab['nk'])]
        for ik in range(tbslab['nk']):
            rdotk = sum(tbslab['kvec'][:,ik] * tbbulk['rvec_all'][:,ir])
            coeff = np.exp(1j * 2 * np.pi * rdotk) / tbbulk['ndegen_all'][ir]
            hr_fourier[rvec_z][ik] += coeff * tbbulk['hr_all'][ir]

    H1 = [np.zeros((nw_match,nw_match), dtype=complex) for ik in range(tbslab['nk'])]
    for ik in range(tbslab['nk']):
        for i in range(nw_match):
            for j in range(nw_match):
                rvec_z = ind_layer_match[j]-ind_layer_match[i]
                if rvec_z in hr_fourier.keys():
                    H1[ik][i,j] += hr_fourier[rvec_z][ik][ind_bulk_match[i],ind_bulk_match[j]]


    H1 = [np.asmatrix(x) for x in H1]
    if init_spin_corrected:
        H2 = [x[np.ix_(ind_slab_match, ind_slab_match)] for x in tbslab['hk_spn']]
    else:
        H2 = [x[np.ix_(ind_slab_match, ind_slab_match)] for x in tbslab['hk']]
    return H1, H2





# def locality_penalty(U, locality_args):
#     penalty = 0.0
#     nk_s = len(U)
#     for ir in locality_args['irmask']:
#         ur_overlap_mat = np.zeros(U[0].shape,dtype=complex)
#         for iks in range(nk_s):
#             rdotk = sum([r * k for r, k in zip(locality_args['rvec_s'][:,ir],
#                         locality_args['ks'][:,iks])])
#             ur_overlap_mat += U[iks] * np.exp(1j*2*np.pi*rdotk)
#         ur_overlap_mat /= nk_s
#         ur_overlap_masked = np.multiply(locality_args['Amask'][ir],
#                                         ur_overlap_mat)
#         penalty += np.sum(np.multiply( np.conj(ur_overlap_mat), ur_overlap_masked)).real
#     return penalty

# def locality_penalty_grad(U, locality_args):
#     Ur_local = {}
#     nk_s = len(U)
#     for ir in locality_args['irmask']:
#         Ur_local[ir] = np.asmatrix(np.zeros(U[0].shape, dtype=complex))
#         for iks in range(nk_s):
#             rdotk = sum([r * k for r, k in zip(locality_args['rvec_s'][:,ir],
#              locality_args['ks'][:,iks])])
#             Ur_local[ir] += U[iks] * np.exp(1j * 2 * np.pi * rdotk)
#         Ur_local[ir] /= nk_s

#     grad = [np.zeros(U[0].shape,dtype=complex) for i in range(len(U))]
#     for iks in range(nk_s):
#         for ir in locality_args['irmask']:
#             rdotk = sum([r * k for r, k in zip(locality_args['rvec_s'][:,ir],
#                         locality_args['ks'][:,iks])])
#             grad[iks] += np.multiply(locality_args['Amask'][ir],
#                                      Ur_local[ir]) \
#                          * np.exp(-1j*2*np.pi*rdotk)
#         grad[iks] /= nk_s
#     return grad



# def run_opt(opt_args):
#     # define variables from input argumentsopt_args
#     H1 = opt_args['H1']
#     H2 = opt_args['H2']
#     U = opt_args['U']
#     iden = opt_args['initial_iden']
#     rvec_s = opt_args['rvec_s']
#     ks = opt_args['ks']
#     nw_s_up = opt_args['nw_s_up']
#     nw_s_dn = opt_args['nw_s_dn']
#     rng_iwb = opt_args['rng_iwb']
#     rng_iws = opt_args['rng_iws']
#     nw_b = opt_args['nw_b']
#     outpath = opt_args['outpath']
#     thr = opt_args['thr']
#     max_iter = opt_args['max_iter']
#     locality_args = opt_args['locality_args']    

#     obj_list = np.empty((0,2))
#     obj_offset = sum([np.linalg.norm(hmat)**2 for hmat in H1]).real
#     obj_offset = np.array([[obj_offset, obj_offset]])
#     obj_list = np.append(obj_list,
#                objective_arr(H1, H2, U, locality_args) + obj_offset, axis=0)
#     norm_Z = thr + 1
#     niter = 0

#     print('nw_s_up = {0:d}, sigma = {1:.0f}'.format(nw_s_up, locality_args['sigma']))
#     f = open(outpath+'output.txt', 'w')
#     f.write('# sigma  niter   norm_Z   error_overlap h_mismatch   penalty\n')
#     f.close()

#     output_args = {}
#     output_args['niter'] = niter
#     output_args['norm_Z'] = norm_Z
#     output_args['obj_list'] = obj_list
#     phase_fix(U, rng_iws, rng_iwb)
#     project_to_unitary(U)
#     opt_args['U'] = U
#     do_postprocess(opt_args, output_args)

#     # main optimization loop
#     mu = 2**-5
#     while (norm_Z > thr):
#         niter += 1
#         if niter > max_iter: break
#         # step 2: calculate gradient D - done inside gradient()
#         # step 3: compute descent direction Z=U*D.H*U-D
#         Z = gradient(H1, H2, U, locality_args)
#         # step 4: evalute norm of Z
#         norm_Z = sum([np.trace(zmat.H*(np.eye(umat.shape[0], dtype=complex)
#                                      -0.5*umat*umat.H)*zmat).real
#                       for umat, zmat in zip(U, Z)])
#         # step 5: update gradient step size mu if too small
#         obj_current = objective(H1, H2, U, locality_args)
#         while True:
#             U_pl_2Z = [umat+2*mu*zmat for umat, zmat in zip(U, Z)]
#             project_to_unitary(U_pl_2Z)
#             dobj_1 = (obj_current
#                       - objective(H1, H2, U_pl_2Z, locality_args))
#             if dobj_1 >= mu * norm_Z:
#                 mu *= 2; print("mu = ", mu)
#             else: break
#             break
#         # step 6: update gradient step size mu if too large
#         while True:
#             U_pl_Z = [umat+mu*zmat for umat, zmat in zip(U, Z)]
#             project_to_unitary(U_pl_Z)
#             dobj_2 = (obj_current
#                       - objective(H1, H2, U_pl_Z, locality_args))
#             if dobj_2 < mu * 0.5 * norm_Z:
#                 mu *= 0.5; print("mu = ", mu)
#             else: break
#             if mu < 1E-10: break
#         if mu < 1E-10: break
#         # setp 7: set U, goto step 2
#         U = U_pl_Z
#         phase_fix(U, rng_iws, rng_iwb)
#         obj_list = np.append(obj_list, 
#             objective_arr(H1, H2, U, locality_args) + obj_offset, axis=0)
#         if niter%10 == 0:
#             print("niter={0}, norm_Z={1:.1E}, dobj_1={2:.1E}, dobj_2={3:.1E}"\
#                   .format(niter, norm_Z, dobj_1, dobj_2))
#         if niter%5 == 0 or niter<100:
#             output_args['niter'] = niter
#             output_args['norm_Z'] = norm_Z
#             output_args['obj_list'] = obj_list
#             opt_args['U'] = U
#             do_postprocess(opt_args, output_args)
    
#     # save U matrices
#     # f = open('figures_GeTe_U_sigma_{0:6.0f}.pkl'.format(sigma), 'wb')
#     # pickle.dump(U, f)
#     # f.close()
    
#     # plot minimization objective and penalty function
#     plt.plot(np.arange(len(obj_list)), obj_list[:,0])
#     plt.plot(np.arange(len(obj_list)), obj_list[:,1]-obj_list[:,0])
#     plt.plot(np.arange(len(obj_list)), obj_list[:,1])
#     plt.legend(['hamiltonian mismatch','penalty','mismatch+penalty'])
#     plt.savefig(outpath+'minimization.png')
#     plt.close()

# def do_postprocess(opt_args, output_args):
#     U = opt_args['U']
#     Uk_overlap = opt_args['Uk_overlap']
#     Ur_overlap = opt_args['Ur_overlap']
#     iden = opt_args['initial_iden']
#     rvec_s = opt_args['rvec_s']
#     ks = opt_args['ks']
#     nw_s_up = opt_args['nw_s_up']
#     nw_s_dn = opt_args['nw_s_dn']
#     nw_b = opt_args['nw_b']
#     outpath = opt_args['outpath']
#     sigma = opt_args['locality_args']['sigma']
#     error_offset = opt_args['locality_args']['error_offset']
#     nk_s = len(U)

#     niter = output_args['niter']
#     norm_Z = output_args['norm_Z']
#     obj_list = output_args['obj_list']

#     Ur = []
#     for ir in range(rvec_s.shape[1]):
#         Ur.append(np.asmatrix(np.zeros(U[0].shape, dtype=complex)))
#         for iks in range(nk_s):
#             rdotk = sum([r * k for r, k in zip(rvec_s[:,ir], ks[:,iks])])
#             Ur[-1] += U[iks] * np.exp(1j * 2 * np.pi * rdotk)
#         Ur[-1] /= nk_s
        
#     ##################
#     #for overlap
#     ##################
#     if niter%100==0 or (niter<100 and niter%10==0) or (niter<10):
#         extent = [-0.5, Ur[0].shape[1]-0.5, Ur[0].shape[0]-0.5, -0.5]
#         plot_irlist = []
#         plot_irlist.append(np.where((rvec_s[0,:] == 0) & (rvec_s[1,:] == 0)
#                             & (rvec_s[2,:] == 0))[0].item())
#         plot_irlist.append(np.where((rvec_s[0,:] == 0) & (rvec_s[1,:] == 1)
#                             & (rvec_s[2,:] == 0))[0].item())
#         plot_irlist.append(np.where((rvec_s[0,:] == 0) & (rvec_s[1,:] == 2)
#                             & (rvec_s[2,:] == 0))[0].item())
#         for ir in plot_irlist:
#             offset = iden.copy()
#             if not (rvec_s[:,ir] == [0,0,0]).all(): offset *= 0
#             fig, axes = plt.subplots(1,3,figsize=(12,6))
#             fig.suptitle(rvec_s[:,ir])
#             vmax = abs(Ur[ir]-offset).max()
#             vmin = -vmax
#             im = axes[0].imshow((Ur[ir]-offset).real, extent=extent, origin='upper',
#                      vmax=vmax, vmin=vmin, cmap='bwr')
#             axes[0].set_title('hamiltonian')
#             axes[1].imshow((Ur_overlap[ir]-offset).real, extent=extent, origin='upper',
#                        vmax=vmax, vmin=vmin, cmap='bwr')
#             axes[1].set_title('overlap')
#             axes[2].imshow((Ur[ir]-Ur_overlap[ir]).real, extent=extent, origin='upper',
#                        vmax=vmax, vmin=vmin, cmap='bwr')
#             axes[2].set_title('ham-overlap')
#             fig.colorbar(im, ax=axes.ravel().tolist(), orientation='horizontal')
#             plt.savefig(outpath+'vmat_ir_{0:d}_niter_{1:05d}.png'.format(ir, niter))
#             plt.close()
            
#     error_overlap = sum([np.linalg.norm(x-y)**2 for x,y in zip(U,Uk_overlap)])
#     error_overlap += error_offset
#     f = open(outpath+'output.txt', 'a')
#     f.write('{0:6.0f} {1:5d} {2:12.4f} {3:10.5f} {4:10.5f} {5:10.5f}\n'.format(sigma,
#             niter, norm_Z, error_overlap, obj_list[-1,0], obj_list[-1,1]-obj_list[-1,0]))
#     f.close()

