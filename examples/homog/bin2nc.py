#!/usr/bin/env python
"""
Usage:
  bin2nc.py plot FILE FIELD
  bin2nc.py report FILE
  bin2nc.py column FILE
  bin2nc.py nc <in> <out>
  bin2nc.py FILE

Arguments:
  FILE                output file
  FIELD               field to plot

"""
import os
import struct
from docopt import docopt
import numpy as np

import matplotlib as mpl


def figsize(scale):
    fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "text.fontsize": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)
import matplotlib.pyplot as plt


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
    """Open bin file"""
    head = read_header(name)
    data = np.memmap(name, dtype=my_dtype(head), offset=head['offset'])

    return head, data


def output2xarray(head, data):
    import xarray as xr

    x = np.arange(head['n']) * head['dx'] / 1000 / 1000
    m = np.arange(1, head['ntrunc'] + 1)
    xticks = [0, 10, 20, 30]
    t = data['time']

    coords = {'x': x, 'm': m, 'time':t}
    xrdict = {}


    for field in data.dtype.fields:
        if field == 'time':
            pass
        elif data[field].shape == (t.shape[0], x.shape[0], m.shape[0]):
            xrdict[field] = (['time', 'x', 'm'], data[field])
        elif data[field].shape == (t.shape[0], x.shape[0]):
            xrdict[field] = (['time', 'x'], data[field])


    return xr.Dataset(xrdict, coords)


def column_report(name):
    import matplotlib.pyplot as plt
    
    xr = output2xarray(*read_output(name))

    fig, axs = plt.subplots(5,1,sharex=True)

    xr.scmt.plot(ax=axs[0])


    for i in range(xr.m.shape[0]):
      xr.u.isel(m=i).plot(ax=axs[1])
      xr.th.isel(m=i).plot(ax=axs[2])


    xr.teb.plot(ax=axs[3])
    xr.q.plot(ax=axs[3])

    xr.hc.plot(ax=axs[4])
    xr.hd.plot(ax=axs[4])
    xr.hs.plot(ax=axs[4])

    axs[0].set_xlim([xr.time.min(), xr.time.max()])

    for ax in axs:
        ax.set_title('')
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    fig.savefig('output.pgf')

    plt.show()


def report(name):
    plt.rc('axes.formatter', limits=(-3, 3), use_mathtext=True)
    hr = 3600.0
    day = hr * 24
    km = 1000

    if os.path.dirname(name) == '':
        name = os.path.join(os.getcwd(), name)
    head, data = read_output(name)

    wd = os.getcwd()
    os.chdir(os.path.dirname(name))

    x = np.arange(head['n']) * head['dx'] / 1000 / 1000
    xticks = [0, 10, 20, 30]
    t = data['time']

    figheight = 9
    width = 7

    # velocity
    fig = plt.figure(1, (width, figheight))
    grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(1, 4),  # creates 2x2 grid of axes
                     cbar_mode='single',
                     aspect=False)

    for i in range(head['ntrunc']):
        cs = grid[i].contourf(x,
                              data['time'],
                              data['u'][..., i],
                              21,
                              cmap='bwr')
        grid.cbar_axes[i].colorbar(cs)
        grid[i].set_title(r'$u_{i}$'.format(i=i + 1))
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

    for i in range(head['ntrunc']):
        cs = grid[i].contourf(x,
                              data['time'],
                              data['th'][..., i],
                              21,
                              cmap='bwr')
        grid.cbar_axes[i].colorbar(cs)
        grid[i].set_title(r'$\theta_{i}$'.format(i=i + 1))
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
        plt.text(.1,
                 .1,
                 field,
                 transform=grid[i].transAxes,
                 bbox=dict(color='w'))
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
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from mpl_toolkits.axes_grid1 import ImageGrid

    args = docopt(__doc__)

    if args['plot']:
        head, data = read_output(args['FILE'])
        from pylab import pcolormesh, show, colorbar, axis

        x = np.arange(head['n']) * head['dx'] / 1000 / 1000
        xticks = [0, 10, 20, 30]
        t = data['time']

        pcolormesh(x, t, data[args['FIELD']])
        axis('tight')
        colorbar()
        show()

    elif args['report']:
        report(args['FILE'])
    elif args['column']:
        column_report(args['FILE'])
    elif args['nc']:
        head, data = read_output(args['<in>'])
        output2xarray(head, data).to_netcdf(args['<out>'])
    else:
        print(read_output(args['FILE'])[0])
