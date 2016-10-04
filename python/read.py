"""Read python simulation output files and create reports

Usage:
    read.py <datadir>
"""
import os
import numpy as np


def read_data(datadir):
    files = [x[:-1]
             for x in open(os.path.join(datadir, "datafiles.txt")).readlines()]
    return np.concatenate([np.load(os.path.join(datadir, fn))['arr_0']
                           for fn in files])

def report_data(data):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.contourf(data['u'][:, 1, :], 12, cmap='bwr')
    plt.colorbar()
    plt.savefig("u.png")

    plt.figure()
    plt.plot(data['time'], np.sqrt(np.mean(data['u']**2, axis=(1, 2))))
    plt.savefig("u_rms.png")

def main():
    from docopt import docopt
    args = docopt(__doc__)

    report_data(read_data(args['<datadir>']))

if __name__ == '__main__':
    main()


