import numpy as np

def find_argument(full_data, keyword, arg_type="text", arr_shape=None):
    '''find a line in data with the structure [keyword = value]
    and return value as output with type arg_type'''
    assert arg_type in ["text", "int", "float", "logical", "int_arr", "float_arr"]

    for line in full_data:
        if line[0] == '#': continue
        if line.startswith(keyword):
            linedata = line.replace(',',' ').replace('=',' ').replace(':',' ').split()
            if "arr" in arg_type:
                value = linedata[1:]
            else:
                value = linedata[-1]

    if arg_type == "text":
        return value
    elif arg_type == "int":
        return int(value)
    elif arg_type == "float":
        return float(value)
    elif arg_type == "logical":
        if value.lower() == "true":
            return True
        elif value.lower() == "false":
            return False
        else:
            raise ValueError
    elif arg_type == "int_arr":
        if arr_shape:
            return np.array(value).reshape(arr_shape).astype(int)
        else:
            return np.array(value).astype(int)
    elif arg_type == "float_arr":
        if arr_shape:
            return np.array(value).reshape(arr_shape).astype(float)
        else:
            return np.array(value).astype(float)


def setup_input(input_filename):
    input_params = {}
    # read input_filename and extract arguments
    with open(input_filename, 'r') as file:
        data = file.readlines()
        material = find_argument(data, 'material')
        path = find_argument(data, 'path')
        
        input_params['seedname_b'] = find_argument(data, 'seedname_b')
        input_params['seedname_s'] = find_argument(data, 'seedname_s')
        input_params['isspinor'] = find_argument(data, 'isspinor', 'logical')

        input_params['iatm_add'] = find_argument(data, 'iatm_add', 'int')

        input_params['nbnd_s'] = find_argument(data, 'nbnd_s', 'int')
        input_params['nw_s'] = find_argument(data, 'nw_s', 'int')
        input_params['nklist_s'] = find_argument(data, 'nklist_s', 'int_arr')
        input_params['nbnd_b'] = find_argument(data, 'nbnd_b', 'int')
        input_params['nw_b'] = find_argument(data, 'nw_b', 'int')
        input_params['nklist_b'] = find_argument(data, 'nklist_b', 'int_arr')

        input_params['iw_head'] = find_argument(data, 'iw_head', 'int')
        input_params['iw_tail'] = find_argument(data, 'iw_tail', 'int')

        input_params['norbital_list'] = list(find_argument(data, 'norbital_list', 'int_arr'))

        input_params['alat_s'] = find_argument(data, 'alat_s', 'float_arr', arr_shape=(3,3)).transpose()
        input_params['alat_b'] = find_argument(data, 'alat_b', 'float_arr', arr_shape=(3,3)).transpose()
        
        return material, path, input_params

