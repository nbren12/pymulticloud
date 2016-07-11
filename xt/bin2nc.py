#!/usr/bin/env python
"""
Usage:
  bin2nc.py plot FILE FIELD
  bin2nc.py report FILE

Arguments:
  FILE                output file
  FIELD               field to plot

"""
import os
import struct
from docopt import docopt
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import ImageGrid



def read_header(name='real.bin'):

    fmt = 'iidddd'
    with open(name, "rb") as f:
        s = struct.unpack_from(fmt, f.read(struct.calcsize(fmt)))

    out = dict(zip(('n', 'ntrunc', 'dx', 'T', 'L', 'abar'), s))
    out['offset'] = struct.calcsize(fmt)

    return out


def my_dtype(head):

    n = head['n']
    ntrunc = head['ntrunc']

    mydt = np.dtype([('time', np.float64),
                     ('u', np.float64, (n, ntrunc)),
                     ('th', np.float64, (n, ntrunc)),
                     ('q', np.float64, (n, )),
                     ('teb', np.float64, (n, )),
                     ('hc', np.float64, (n, )),
                     ('hd', np.float64, (n, )),
                     ('hs', np.float64, (n, )),
                     ('fcls', np.float64, (n, )),
                     ('fdls', np.float64, (n, )),
                     ('fsls', np.float64, (n, )),
                     ('scmt', np.int32, (n, ))])

    return mydt


def read_output(name):
    head = read_header(name)
    data = np.memmap(name, dtype=my_dtype(head), offset=head['offset'])

    return head, data


def report(name):
    plt.rc('axes.formatter', limits=(-3,3), use_mathtext=True)
    hr = 3600.0
    day = hr * 24
    km = 1000

    head, data = read_output(name)


    wd = os.getcwd()
    os.chdir(os.path.dirname(name))

    x = np.arange(head['n']) * head['dx']/ 1000/1000
    xticks = [0, 10, 20, 30]
    t = data['time']

    figheight= 9
    width = 7

    # velocity
    fig = plt.figure(1, (width, figheight))
    grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(1, 4),  # creates 2x2 grid of axes
                     cbar_mode='single',
                     aspect=False)


    for i in range(4):
        cs = grid[i].contourf(x, data['time'], data['u'][...,i], 21, cmap='bwr')
        grid.cbar_axes[i].colorbar(cs)
        grid[i].set_title(r'$u_{i}$'.format(i=i+1))
        grid[i].set_xticks(xticks)

    plt.tight_layout()
    plt.subplots_adjust(right=.94)
    fig.savefig("velocity.png")


    # temperature
    fig = plt.figure(2, (width, figheight))
    grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(1, 4),  # creates 2x2 grid of axes
                     cbar_mode='single',
                     aspect=False)

    for i in range(4):
        cs = grid[i].contourf(x, data['time'], data['th'][...,i], 21, cmap='bwr')
        grid.cbar_axes[i].colorbar(cs)
        grid[i].set_title(r'$\theta_{i}$'.format(i=i+1))
        grid[i].set_xticks(xticks)

    plt.tight_layout()
    plt.subplots_adjust(right=.91, top=.96)
    fig.savefig("temp.png")

    # thermo
    fig = plt.figure(3, (width, figheight))
    grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(1, 5),  # creates 2x2 grid of axes
                     cbar_mode='each',
                     cbar_location='top',
                     aspect=False)

    for i, field in enumerate(['teb', 'q', 'hc', 'hd', 'hs']):
        if field in ['hc', 'hd', 'hs']:
            cmap = 'YlGnBu_r'
        else:
            cmap = 'bwr'
        cs = grid[i].contourf(x, data['time'], data[field], 21, cmap=cmap)
        cb = grid.cbar_axes[i].colorbar(cs)
        cb.locator = plt.MaxNLocator(4)
        plt.text(.1, .1, field, transform=grid[i].transAxes, bbox=dict(color='w'))
        grid[i].set_xticks(xticks)

    plt.tight_layout()
    plt.subplots_adjust(top=.96)
    fig.savefig("thermo.png")

    plt.figure(figsize=(2, figheight))

    plt.pcolormesh(x, data['time'], data['scmt'], cmap=cmap)
    plt.gca().set_xticks(xticks)
    plt.axis('tight')
    plt.savefig("scmt.png")

    os.chdir(wd)



if __name__ == '__main__':

    args = docopt(__doc__)

    if args['plot']:
        head, data = read_output(args['FILE'])
        from pylab import pcolormesh, show, colorbar, axis

        x = np.arange(head['n']) * head['dx']/ 1000/1000
        xticks = [0, 10, 20, 30]
        t = data['time']

        pcolormesh(x, t, data[args['FIELD']])
        axis('tight')
        colorbar()
        show()

    if args['report']:
        report(args['FILE'])

