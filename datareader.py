import h5py

def hdf5_to_dict(group):
    dict = {}
    for key, item in group.items():
        if isinstance(item, h5py.Group):
            dict[key] = hdf5_to_dict(item)
        else:
            dict[key] = item[()]
    return dict

def extract(path):
    
    if path!="":
        ex_data = h5py.File(path+'Fields0.h5', 'r')
        return hdf5_to_dict(ex_data)
    else:
        ex_data = h5py.File('Fields0.h5', 'r')
        return hdf5_to_dict(ex_data)

data_dict = extract("../azim_4/") 