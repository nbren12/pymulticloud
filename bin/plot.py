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

from python.read import read_data, read_diags
from python.cmt import calc_du

data = read_data("./data")
diags = read_diags("./diags.pkl")


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

    next(pt)
    plt.plot(d['ncmt'])
    plt.title("ncmt")
    plt.legend(range(3))

    next(pt)
    plt.title("kcmt")
    plt.plot(d['kcmt'])

    next(pt)
    plt.title("kcmt * ncmt")
    plt.plot(d['kcmt'] * d['ncmt'])
    plt.legend(range(3))

    next(pt)
    plt.title("rms")
    plt.plot(d['rms'])

    next(pt)
    rms = np.ma.array(d['rms'])
    n = np.sum(d['ncmt'][1, :])
    fcmt = np.sum(d['kcmt'] * d['ncmt'] / n, axis=1)
    keff = fcmt / rms
    plt.plot(keff)
    plt.title(r"$k_{{eff}}^{{-1}}$ = 1/{:.2f} = {:.2f} day".format(-keff.mean(
    ), -1 / keff.mean() * (8.3 / 24)))

    plt.tight_layout()

with PdfPages("report.pdf") as pdf:

    def clim_plot(data):
        # Climatology plots
        variables = list(data.dtype.fields.keys())
        variables.remove("time")

        if os.path.exists("./averages.npz"):
            avg = np.load("./averages.npz")['avg']
            clim = {key: data[key] for key in variables}
        else:
            clim = {key: np.mean(data[key], axis=0) for key in variables}

        plots = OrderedDict()
        for k, fields in {'fracs': ['fc', 'fd', 'fs'],
                          'rates': ['hc', 'hd', 'hs'],
                          'moist': ['teb', 'q', 'lmd', 'tebst']}.items():
            plots[k] = {f: clim.pop(f) for f in fields}

        # velocity and temperature plots
        for i in range(1, 3):
            plots.setdefault('u', {})['u{}'.format(i)] = clim['u'][i, :]
            plots.setdefault('t', {})['t{}'.format(i)] = clim['t'][i, :]
        clim.pop('t')
        clim.pop('u')

        # any extra plots
        plots['misc'] = clim

        for k in plotiter(plots, w=6, ncol=1, aspect=.3):
            for field in plots[k]:
                data = plots[k][field]
                plt.plot(np.squeeze(data), label=field)
            plt.legend()
        plt.tight_layout()

    def hov_plot(data):

        plots = OrderedDict()
        istart = 400
        plots['u1'] = (data['u'][istart:, 1, :], ), dict(cmap='bwr')
        plots['hd'] = (data['hd'][istart:, ...], ), dict(norm=LogP1(3),
                                                         cmap='YlGnBu_r')
        plots['q'] = (data['q'][istart:, ...], ), dict(cmap='YlGnBu')
        plots['tem-teb'] = (data['lmd'][istart:, ...], ), dict(cmap='YlGnBu')

        if 'scmt' in data.dtype.fields:
            plots['scmt'] = (data['scmt'][istart:, ...], ), dict(
                cmap='YlGnBu_r')

        for k in plotiter(plots,
                          w=4,
                          ncol=3,
                          aspect=2,
                          label_dict=dict(loc=(0.0, 1.03))):
            arg, kw = plots[k]
            plt.pcolormesh(*arg, rasterized=True, **kw)
            plt.colorbar()
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

    if 'kcmt' in diags:
        kcmt_plots(diags)
        pdf.savefig()
