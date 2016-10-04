"""Implementation of multicloud model in python

This code uses the forcings and gillespie algorithm from the fortran code
"""
import itertools
import sys
import os
import logging
import uuid

import numpy as np

from .two_mode_swe import f as f2m
from .tadmor_1d import periodic_bc, central_scheme
from .timestepping import steps

# TODO add the multicloud code to this project
from fortran import multicloud as mc

L = 3

variables = ['q', 'teb', 'hs', 'tebst', 'fc', 'fd', 'fs']
variable_idxs = dict(zip(variables, itertools.count(2 * L)))
variable_idxs['u'] = slice(0, 2 * L, 2)
variable_idxs['t'] = slice(1, 2 * L, 2)


def f(q, alpha_tld=0.1, lmd_tld=0.8, q_tld=0.9):
    """Flux function for multicloud model"""

    u = q[variable_idxs['u'], ...]
    T = q[variable_idxs['t'], ...]
    moist = q[variable_idxs['q'], ...]

    fq = np.empty_like(q)
    f2m(q, fq=fq)

    # moisture equation
    fq[variable_idxs['q'], ...] = q_tld * (u[1] + lmd_tld * u[2]) \
                                  + moist * (u[1] + alpha_tld * u[2])

    return fq


def onestep(soln, time, dt, dx):

    # hyperbolic terms
    periodic_bc(soln)
    soln[:2 * L + 1, ...] = central_scheme(f, soln[:2 * L + 1, ...], dx, dt)

    # multicloud model step
    mc.multicloud_rhs(soln[variable_idxs['fc']], soln[variable_idxs['fd']],
                      soln[variable_idxs['fs']], soln[variable_idxs['u']][1],
                      soln[variable_idxs['u']][2], soln[variable_idxs['t']][1],
                      soln[variable_idxs['t']][2], soln[variable_idxs['teb']],
                      soln[variable_idxs['q']], soln[variable_idxs['hs']], dt,
                      dx, time, soln[variable_idxs['tebst']])

    return soln


def init_mc(n=1000, dx=50 / 1500):
    import numpy as np
    neq = 2 * L + len(variables)
    soln = np.zeros((neq, n))

    fceq, fdeq, fseq = mc.equilibrium_fractions()

    soln[variable_idxs['fc']] = fceq
    soln[variable_idxs['fd']] = fdeq
    soln[variable_idxs['fs']] = fseq

    # initialize temperature field with small random perturbation
    soln[variable_idxs['t'], ...][0] = np.random.randn(n) * .01

    return soln, dx


def init_mc_from_file(fn):
    neq = 2 * L + len(variables)
    icdata = np.load(fn)['arr_0'][-1]
    n = icdata['teb'].shape[0] + 4
    soln = np.zeros((neq, n))

    for key in variable_idxs:
        soln[variable_idxs[key]][..., 2:-2] = icdata[key]

    return soln


def unghosted(q):
    return q[:, 2:-2]


def record_array_soln(soln, t):
    """Create record array for solution"""

    soln = unghosted(soln)
    n = soln.shape[-1]

    variables_dtype = [(k, np.float64, (n,)) for k in variables] \
                      + [(k, np.float64, (L, n)) for k in ['u', 't']]
    mydt = np.dtype([('time', np.float64)] + variables_dtype)
    arr = np.zeros(1, dtype=mydt)

    for k, v in variable_idxs.items():
        arr[k] = soln[v]

    arr['time'] = t

    return arr


def save_restart_file(name, soln, t, dx):
    import pickle
    logging.info("Saving restart file at t={0}".format(t))
    with open(name, "wb") as f:
        pickle.dump((soln, t, dx), f)


def load_restart_file(name):
    import pickle
    with open(name, "rb") as f:
        arr = pickle.load(f)

    return arr


def main():
    """Runs multicloud model

    TODO This file is too complicated needs to be refactored, and the IO needs
    to be rethought
    """
    logging.basicConfig(level=logging.INFO)

    soln, dx = init_mc()
    t_start = 0.0
    t_run = 100

    if os.path.exists('restart.pkl'):
        soln, t_start, dx = load_restart_file("restart.pkl")
        logging.info("Loading restart file at t={0}".format(t_start))
    elif os.path.exists('ic.npz'):
        soln = init_mc_from_file("ic.npz")

    dt = dx * .1
    t_end = t_start + t_run

    dt_out = 1.0
    t_out = t_start + dt_out
    i_out = 1
    i_file = 0

    arr = record_array_soln(soln, 0.0)

    nbuf = 10
    output = np.zeros(nbuf, dtype=arr.dtype)
    output[0] = arr

    datadir = "data"

    if not os.path.isdir(datadir):
        os.mkdir(datadir)
    f = open(os.path.join(datadir, "datafiles.txt"), "w+")

    for t, soln in steps(onestep, soln, dt, (t_start, t_end), dx):

        if t > t_out:
            logging.info("Storing output data at t={0}".format(t))
            output[i_out] = record_array_soln(soln, t)
            t_out += dt_out
            i_out += 1

            if i_out >= nbuf:
                cur_file_name = str(uuid.uuid1()) + ".npz"
                logging.info("Saving data to file `{1}` at t={0}".format(
                    t, cur_file_name))
                np.savez(os.path.join(datadir, cur_file_name), output)
                f.write(cur_file_name)
                f.write("\n")
                i_file += 1
                i_out = 0

    f.close()
    save_restart_file("restart_" + str(uuid.uuid1()) + ".pkl", soln, t, dx)
    np.savez("soln%05d.npz" % i_file, output[:i_out])


def test_record_array_soln():

    soln, dx, dt = init_mc()
    arr = record_array_soln(soln, 0.0)


if __name__ == '__main__':
    main()
