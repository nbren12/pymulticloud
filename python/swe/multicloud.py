"""Implementation of multicloud model in python

This code uses the forcings and gillespie algorithm from the fortran code
"""
import itertools
import sys
import os
import logging
import uuid

import numpy as np

# this import needs to happen before the others for some reason. This is
# probably a conflict with numba.
import fortran.multicloud as mc

from .two_mode_swe import f as f2m
from .tadmor_1d import periodic_bc, central_scheme
from .timestepping import steps

logger = logging.getLogger(__file__)


class MulticloudModel(object):
    L = 3
    variables = ['q', 'teb', 'hs', 'tebst', 'fc', 'fd', 'fs', 'hc', 'hd']

    @property
    def variable_idxs(self):
        variable_idxs = dict(zip(self.variables, itertools.count(2 * self.L)))
        variable_idxs['u'] = slice(0, 2 * self.L, 2)
        variable_idxs['t'] = slice(1, 2 * self.L, 2)

        return variable_idxs

    def _f(self, q, alpha_tld=0.1, lmd_tld=0.8, q_tld=0.9, nonlinear=0.0):
        """Flux function for multicloud model"""
        variable_idxs = self.variable_idxs

        u = q[variable_idxs['u'], ...]
        T = q[variable_idxs['t'], ...]
        moist = q[variable_idxs['q'], ...]

        fq = np.empty_like(q)
        f2m(q, fq=fq, nonlin=nonlinear)

        # moisture equation
        fq[variable_idxs['q'], ...] = q_tld * (u[1] + lmd_tld * u[2]) \
                                    + moist * (u[1] + alpha_tld * u[2])

        return fq

    def onestep(self, soln, time, dt, dx, nonlinear=1.0):
        """Perform a single time step of the multicloud model"""
        from functools import partial
        variable_idxs = self.variable_idxs

        # hyperbolic terms
        periodic_bc(soln)
        f_partial = partial(self._f, nonlinear=nonlinear)
        soln[:2 * self.L + 1, ...] = central_scheme(
            f_partial, soln[:2 * self.L + 1, ...], dx, dt)

        # multicloud model step
        mc.multicloud_rhs(
            soln[variable_idxs['fc']], soln[variable_idxs['fd']],
            soln[variable_idxs['fs']], soln[variable_idxs['u']][1],
            soln[variable_idxs['u']][2], soln[variable_idxs['t']][1],
            soln[variable_idxs['t']][2], soln[variable_idxs['teb']],
            soln[variable_idxs['q']], soln[variable_idxs['hs']], dt, dx, time,
            soln[variable_idxs['tebst']], soln[variable_idxs['hc']],
            soln[variable_idxs['hd']])

        return soln

    @property
    def neq(self):
        return 2 * self.L + len(self.variables)

    def init_mc(self, n=1000, dx=40 / 1500, asst=0.5, lsst=10000 / 1500):
        variable_idxs = self.variable_idxs

        soln = np.zeros((self.neq, n))

        fceq, fdeq, fseq = mc.equilibrium_fractions()

        soln[variable_idxs['fc']] = fceq
        soln[variable_idxs['fd']] = fdeq
        soln[variable_idxs['fs']] = fseq

        # initialize temperature field with small random perturbation
        soln[variable_idxs['t'], ...][0] = np.random.randn(n) * .01

        # initialize teb
        x = np.arange(n) * dx
        domain_size = n * dx

        tebst = asst * np.cos(2 * np.pi * (x - domain_size / 2) / lsst) * (
            x > domain_size / 2 - lsst / 2) * (x < domain_size / 2 + lsst / 2)
        tebst[x >= domain_size / 2 + lsst / 2] = -asst
        tebst[x <= domain_size / 2 - lsst / 2] = -asst

        soln[variable_idxs['tebst']] = tebst

        return soln, dx

    def init_mc_from_file(self, fn):

        variable_idxs = self.variable_idxs
        neq = self.neq
        icdata = np.load(fn)['arr_0'][-1]
        n = icdata['teb'].shape[0] + 4
        soln = np.zeros((neq, n))

        for key in variable_idxs:
            soln[variable_idxs[key]][..., 2:-2] = icdata[key]

        return soln

    def record_array_soln(self, soln, t):
        """Create record array for solution"""
        variables = self.variables
        variable_idxs = self.variable_idxs

        soln = unghosted(soln)
        n = soln.shape[-1]

        variables_dtype = [(k, np.float64, (n,)) for k in variables] \
                        + [(k, np.float64, (self.L, n)) for k in ['u', 't']]
        mydt = np.dtype([('time', np.float64)] + variables_dtype)
        arr = np.zeros(1, dtype=mydt)

        for k, v in variable_idxs.items():
            arr[k] = soln[v]

        arr['time'] = t

        return arr


def unghosted(q):
    return q[:, 2:-2]


def save_restart_file(name, soln, t, dx):
    import pickle
    logger.info("Saving restart file `{file}` at t={0}".format(t, file=name))
    with open(name, "wb") as f:
        pickle.dump((soln, t, dx), f)


def load_restart_file(name):
    import pickle
    with open(name, "rb") as f:
        arr = pickle.load(f)

    return arr


def main(run_duration=100, dt_out=1.0, solver=None):
    """Runs multicloud model

    TODO This file is too complicated needs to be refactored, and the IO needs
    to be rethought

    Parameters
    ----------
    run_duration: float
        length of simulation from start_time
    dt_out: float
        output interval
    """
    t_start = 0.0

    logger.info("Starting run with duration={0}".format(run_duration))

    if solver is None:
        solver = MulticloudModel()

    soln, dx = solver.init_mc()

    if os.path.exists('restart.pkl'):
        soln, t_start, dx = load_restart_file("restart.pkl")
        logger.info("Loading restart file at t={0}".format(t_start))
    elif os.path.exists('ic.npz'):
        soln = solver.init_mc_from_file("ic.npz")

    dt = dx * .1
    t_end = t_start + run_duration

    t_out = t_start + dt_out

    # allocate output buffer
    nbuf = 10
    arr = solver.record_array_soln(soln, t_start)
    output = np.zeros(nbuf, dtype=arr.dtype)

    # include initial data
    output[0] = arr
    i_out = 1

    datadir = "data"

    if not os.path.isdir(datadir):
        os.mkdir(datadir)

    for t, soln in steps(solver.onestep, soln, dt, (t_start, t_end), dx):

        if t > t_out:
            logger.info("Storing output data at t={0}".format(t))
            output[i_out] = solver.record_array_soln(soln, t)
            t_out += dt_out
            i_out += 1

            if i_out >= nbuf:
                dump_output_file(t, output, datadir)
                i_out = 0

    dump_output_file(t, output[:i_out], datadir)
    save_restart_file("restart_" + str(uuid.uuid1()) + ".pkl", soln, t, dx)


def dump_output_file(t, output, datadir):

    cur_file_name = str(uuid.uuid1()) + ".npz"

    logger.info("Saving data to file `{1}` at t={0}".format(t, cur_file_name))
    np.savez(os.path.join(datadir, cur_file_name), output)

    with open(os.path.join(datadir, "datafiles.txt"), "a") as f:
        f.write(cur_file_name)
        f.write("\n")


def test_record_array_soln():

    soln, dx, dt = init_mc()
    arr = record_array_soln(soln, 0.0)


if __name__ == '__main__':
    main()
