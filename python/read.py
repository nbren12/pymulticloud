"""Read python simulation output files and create reports

Usage:
    read.py <datadir>
"""
import os
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
