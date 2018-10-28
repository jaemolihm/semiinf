import numpy as np

def setup_input():
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
    return material, path, input_params