#!/usr/bin/env python
""" Plotting script for python based multicloud model runs
"""
import os
import sys

import numpy as np
from gnl.plots import *



root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, root_dir)

from python.read import read_data, read_diags
from python.cmt import calc_du



data = read_data("./data")
diags = read_diags("./diags.pkl")

def clim_plot(data):
    variables = list(data.dtype.fields.keys())

    variables.remove("u")
    variables.remove("t")
    variables.remove("time")
    clim = {key: np.mean(data[key], axis=0) for key in variables}
    for k, v in plotiter(clim.items(), w=3, aspect=.8):
        plt.plot(v)
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
    plt.plot(d['kcmt']*d['ncmt'])
    plt.legend(range(3))

    next(pt)
    plt.title("rms")
    plt.plot(d['rms'])

    next(pt)
    rms = np.ma.array(d['rms'])
    n = np.sum(d['ncmt'][1,:])
    fcmt = np.sum(d['kcmt']*d['ncmt']/n, axis=1)
    keff=  fcmt/rms
    plt.plot(keff)
    plt.title(r"$k_{{eff}}^{{-1}}$ = 1/{:.2f} = {:.2f} day".format(-keff.mean(), -1/keff.mean() * (8.3/24)))

    plt.tight_layout()

with PdfPages("report.pdf") as pdf:

    clim_plot(data)
    pdf.savefig()
    plt.close()

    aspect = data['u'].shape[0] / 100 * .3
    # cmt plot
    if 'scmt' in diags:
        for i in plotiter(range(1,3), w=3, aspect=aspect):
            plt.pcolormesh(data['scmt']==i, cmap='Greys', rasterized=True)
            plt.title("CMT={0}".format(i))
        pdf.savefig()
        plt.close()

    # u plot
    for i in plotiter(range(1,3), w=3, aspect=aspect):
        plt.pcolormesh(data['u'][:,i,:], cmap='bwr', rasterized=True)
        plt.title(r"$u_{}$".format(i))
    pdf.savefig()
    plt.close()

    plt.figure()
    umean = data['u'][:,:,2:-2].mean(axis=0)
    plt.plot(umean.T)
    plt.legend([0,1,2])
    pdf.savefig()

    plt.figure()
    dul, duhi = calc_du(umean)
    plt.plot(dul, label='dulow')
    plt.plot(duhi, label='duhi')
    plt.legend()
    pdf.savefig()


    plt.figure()
    plt.plot(np.sqrt(data['u'][:,:,2:-2]**2).mean(axis=0).T)
    plt.legend([0,1,2])
    plt.grid()
    pdf.savefig()


    plt.figure(figsize=(4,4*aspect))
    plt.pcolormesh(data['hd'], rasterized=True, norm=LogP1(3), cmap='YlGnBu_r')
    plt.title(r'$H_d$')
    plt.colorbar()
    pdf.savefig()

    if 'kcmt' in diags:
        kcmt_plots(diags)
        pdf.savefig()

