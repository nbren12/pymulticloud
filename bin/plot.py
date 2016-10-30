#!/usr/bin/env python
""" Plotting script for python based multicloud model runs
"""
import os
import sys

from collections import OrderedDict
import numpy as np
from gnl.plots import *

root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, root_dir)

from stochmc.read import read_data, read_diags, load_xarray
from stochmc.cmt import calc_du

data  = load_xarray("data/", "diags.pkl")


def clim_plot_mean(data):
    variables = list(data.dtype.fields.keys())

    variables.remove("u")
    variables.remove("t")
    variables.remove("time")
    clim = {key: np.mean(data[key], axis=0) for key in variables}
    for k, v in plotiter(clim.items(), w=3, aspect=.8):
        plt.plot(v)
        plt.title(k)


def clim_plot(data):
    variables = list(data.dtype.fields.keys())

    variables.remove("u")
    variables.remove("t")
    variables.remove("time")
    clim = {key: data[key] for key in variables}
    for k, v in plotiter(clim.items(), w=3, aspect=.8):
        plt.plot(np.squeeze(v))
        plt.title(k)


def kcmt_plots(d):
    pt = plotiter(range(5), ncol=1, w=6, aspect=.3)

    kcmt = np.vstack(d['kcmt{}'.format(i)] for i in range(3)).T
    ncmt = np.vstack(d['ncmt{}'.format(i)] for i in range(3)).T

    next(pt)
    plt.plot(ncmt)
    plt.title("ncmt")
    plt.legend(range(3))

    next(pt)
    plt.plot(kcmt)
    plt.title("kcmt")
    plt.legend(range(3))

    next(pt)
    plt.title("kcmt * ncmt")
    plt.plot(kcmt * ncmt)
    plt.legend(range(3))

    next(pt)
    plt.title("rms")
    plt.plot(d['rms'])

    next(pt)
    rms = np.ma.array(d['rms'])
    n = np.sum(ncmt[1, :])
    fcmt = np.sum(kcmt * ncmt / n, axis=1)
    keff = fcmt / rms
    plt.plot(keff)
    plt.title(r"$k_{{eff}}^{{-1}}$ = 1/{:.2f} = {:.2f} day".format(-keff.mean(
    ), -1 / keff.mean() * (8.3 / 24)))

    plt.tight_layout()

with PdfPages("report.pdf") as pdf:

    def clim_plot(data):
        # Climatology plots
        clim = data.mean('time')
        variables = set([k for k in clim.data_vars if "i" in data[k].coords])

        plots = OrderedDict()
        for k, fields in {'fracs': ['fc', 'fd', 'fs'],
                          'rates': ['hc', 'hd', 'hs'],
                          'moist': ['teb', 'q', 'lmd', 'tebst']}.items():
            plots[k] = {}
            for f in fields:
                plots[k][f] = clim[f]
                variables.remove(f)

        
        # velocity and temperature plots
        plots['u'] = {}
        plots['t'] = {}
        for i in range(0, 3):
            v = 'u{}'.format(i)
            t = 't{}'.format(i)
            if i > 0:
                plots['u'][v] = clim[v]
                plots['t'][t] = clim[t]
            variables.remove(v)
            variables.remove(t)

        # any extra plots
        if len(variables) > 0:
            plots['misc'] = {}
            for k in variables:
                plots['misc'][k] = clim[k]

        for k in plotiter(plots, w=6, ncol=1, aspect=.3):
            for field in plots[k]:
                data = plots[k][field]
                data.plot()
            plt.legend()
        plt.tight_layout()

    def hov_plot(data):
        from scipy.ndimage import gaussian_filter

        plots = OrderedDict()

        data = data.isel(time=slice(-800, None))


        plots['u1'] = data['u1'], dict(cmap='bwr')


        hdblur = data['hd'].copy()
        hdblur.values = gaussian_filter(hdblur.values, 1.2)
        plots['hd'] = hdblur, dict(cmap='hot', norm=LogP1(data['hd']), vmin=None, vmax=None)
        plots['q'] = data['q'], dict(cmap='YlGnBu')
        plots['tem-teb'] = data['lmd'], dict(cmap='YlGnBu')

        if 'scmt' in data:
            plots['scmt'] = data['scmt'], dict(cmap='bwr', vmin=0, vmax=2)

        for k in plotiter(plots,
                          w=4,
                          ncol=3,
                          aspect=2,
                          label_dict=dict(loc=(0.0, 1.03))):
            arg, kw1 = plots[k]

            vmin, vmax = np.percentile(arg, [1,99])

            kw = {}
            kw['vmin'] = vmin
            kw['vmax'] = vmax
            kw.update(kw1)

            im = arg.plot(rasterized=True, **kw)
            plt.title(k)
            plt.axis('tight')

    clim_plot(data)
    pdf.savefig()
    plt.close()

    hov_plot(data)
    pdf.savefig()
    plt.close()
    # # u plot
    # for i in plotiter(range(1,3), w=3, aspect=aspect):
    #     plt.pcolormesh(data['u'][:,i,:], cmap='bwr', rasterized=True)
    #     plt.title(r"$u_{}$".format(i))
    # pdf.savefig()
    # plt.close()

    # plt.figure()
    # umean = data['u'][:,:,2:-2].mean(axis=0)
    # plt.plot(umean.T)
    # plt.legend([0,1,2])
    # pdf.savefig()

    # plt.figure()
    # dul, duhi = calc_du(umean)
    # plt.plot(dul, label='dulow')
    # plt.plot(duhi, label='duhi')
    # plt.legend()
    # pdf.savefig()

    # plt.figure()
    # plt.plot(np.sqrt(data['u'][:,:,2:-2]**2).mean(axis=0).T)
    # plt.legend([0,1,2])
    # plt.grid()
    # pdf.savefig()

    # plt.figure(figsize=(4,4*aspect))
    # plt.pcolormesh(data['hd'], rasterized=True, norm=LogP1(3), cmap='YlGnBu_r')
    # plt.title(r'$H_d$')
    # plt.colorbar()
    # pdf.savefig()

    if 'kcmt1' in data:
        kcmt_plots(data)
        pdf.savefig()
