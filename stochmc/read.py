"""Read python simulation output files and create reports

Usage:
    read.py <datadir>
"""
import os
import itertools
import numpy as np


def data_iter(datadir):
    files = [os.path.join(datadir, x[:-1])
             for x in open(os.path.join(datadir, "datafiles.txt")).readlines()]
    return (np.load(f)['arr_0'] for f in files)

def read_data(datadir):
    files = [x[:-1]
             for x in open(os.path.join(datadir, "datafiles.txt")).readlines()]
    return np.concatenate([np.load(os.path.join(datadir, fn))['arr_0']
                           for fn in files])
def read_diags(name):
    import pickle
    with open(name, "rb") as f:
        diags = []
        while True:
            try:
                diags.append(pickle.load(f))
            except EOFError:
                break

    out = {}
    for key in diags[0]:
        out[key] = np.vstack([d[key] for d in diags])
        out[key] = np.squeeze(out[key])

    return out


def separate_diags(diags):

    out_dict = {}
    for key in diags:
        arr = np.squeeze(diags[key])

        if arr.ndim > 1:
            axes = arr.shape[1:]

            for inds in itertools.product(*[range(n) for n in axes]):
                key_expanded = key +'x'.join(map(str, inds))

                inds = [slice(None)] + list(inds)
                out_dict[key_expanded] = diags[key][inds]
        else:
            out_dict[key] = diags[key]

    return out_dict

def create_xarray(data, diags):
    import xarray as xr

    scalar_variables = list(data.dtype.fields.keys())

    scalar_variables.remove("u")
    scalar_variables.remove("t")

    scalar_variables.remove("time")

    data_vars = {}
    for key in scalar_variables:
        data_vars[key] = (['time', 'i'], data[key])

    i = np.arange(data['u'].shape[-1])
    coords = {'time': data['time'],
              'i': i}

    # non scalar vars
    for key in ['u', 't']:
        nm = data[key].shape[1]
        for m in range(nm):
            data_vars[key + str(m)] = (['time', 'i'], data[key][:, m,:])

    # diagnostic variables
    diags = separate_diags(diags)
    for key in diags:
        if key != 'time':
            data_vars[key] = (['time_d',], diags[key])

    coords['time_d'] = diags['time']




    return xr.Dataset(data_vars, coords)


def load_xarray(datadir="data/", diagname="diags.pkl"):
    data = read_data(datadir)
    diags = read_diags(diagname)

    return create_xarray(data, diags)

def main(datadir, diagname, output_name):
    data = read_data(datadir)
    diags = read_diags(diagname)

    create_xarray(data, diags).to_netcdf(output_name)

    return 0

if __name__ == '__main__':
    import sys
    main(*sys.argv[1:])
