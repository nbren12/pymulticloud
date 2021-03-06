"""Implementation of multicloud model in python

This code uses the forcings and gillespie algorithm from the fortran code
"""
import itertools
import sys
import os
import pickle
import logging
import uuid
import copy
from collections import deque
import numpy as np
from math import sqrt

from gnl.pdes.tadmor.tadmor import Tadmor1DBase
from gnl.pdes.timestepping import steps

# this import needs to happen before the others for some reason. This is
# probably a conflict with numba.
# from ..wrapper import multicloud as mc
from fortran import mc
from ..io import NetCDF4Writer

logger = logging.getLogger(__file__)

class Soln(object):
    """Object for storing the state of the system
    """
    L = 3
    variables = ['q', 'teb', 'hs', 'tebst', 'fc', 'fd', 'fs', 'hc', 'hd', 'lmd']

    def __init__(self, n, extra_vars=[]):
        "docstring"
        self.variables += list(extra_vars)
        self._data = np.zeros((self.neq, n))

    @classmethod
    def from_arr(cls, arr, **kwargs):
        soln = cls(arr.shape[-1], **kwargs)
        soln._data = arr

        return soln

    def __getitem__(self, name):
        try:
            return self._data[self.variable_idxs[name]]
        except:
            return self._data[name]

    def __setitem__(self, name, val):
        try:
            self._data[self.variable_idxs[name]] = val
        except:
            self._data[name] = val

    # def __getattr__(self, name):
    #     return getattr(self._data, name)

    def comm(self):
        uc = self._data
        uc[:, :2] = uc[:, -4:-2]
        uc[:, -2:] = uc[:, 2:4]


    @property
    def vars(self):
        """Iterator of variable dimension information

        This is useful for output to netcdf.
        """

        yield 'u', ['m', 'x'], self['u'][:,2:-2]
        yield 't', ['m', 'x'], self['t'][:,2:-2]

        for v in self.variables:
            yield v, ['x'], self[v][2:-2]

        return

    @property
    def q(self):
        return self._data[:2*self.L+1,:]

    @q.setter
    def q(self, val):
        self._data[:2*self.L+1,:] = val

    @property
    def neq(self):
        return 2 * self.L + len(self.variables)

    @property
    def variable_idxs(self):
        variable_idxs = dict(zip(self.variables, itertools.count(2 * self.L)))
        variable_idxs['u'] = slice(0, 2 * self.L, 2)
        variable_idxs['t'] = slice(1, 2 * self.L, 2)

        return variable_idxs

    def record_array_soln(self,t):
        """Create record array for solution"""
        variables = self.variables
        variable_idxs = self.variable_idxs

        n = self._data.shape[-1]

        variables_dtype = [(k, np.float64, (n,)) for k in variables] \
                        + [(k, np.float64, (self.L, n)) for k in ['u', 't']]
        mydt = np.dtype([('time', np.float64)] + variables_dtype)
        arr = np.zeros(1, dtype=mydt)

        for k, v in variable_idxs.items():
            arr[k] = self[v]

        arr['time'] = t

        return arr

def f(q, fq=None, alpha_tld=0.1, lmd_tld=0.8, q_tld=0.9, L=3):
    """Conservative flux according to stechmann and majda"""
    u = q[:2*L:2]
    T = q[1:2*L:2]
    moist = q[L*2]

    if fq is None:
        fq = np.empty_like(q)

    fq[0] = 0.0
    fq[1] = 0.0
    fq[2] = - T[1]
    fq[3] = - u[1]
    fq[4] = -T[2]
    fq[5] = -u[2] / 4
    fq[6] = q_tld * (u[1] + lmd_tld * u[2]) \
                                  + moist * (u[1] + alpha_tld * u[2])

    return fq

def f_nonlinear(q, fq=None, alpha_tld=0.1, lmd_tld=0.8, q_tld=0.9, L=3):
    """Conservative flux according to stechmann and majda"""
    u = q[:2*L:2]
    T = q[1:2*L:2]
    moist = q[L*2]

    if fq is None:
        fq = np.empty_like(q)

    fq[0] = 0.0
    fq[1] = 0.0
    fq[2] = - T[1] +3.0/sqrt(2) * u[1] * u[2]
    fq[3] = - u[1] + sqrt(2) *u[1] * T[2] - u[2] * T[1] / sqrt(2)
    fq[4] = -T[2]
    fq[5] = -u[2] / 4
    fq[6] = q_tld * (u[1] + lmd_tld * u[2]) \
                                  + moist * (u[1] + alpha_tld * u[2])

    return fq


def nonlinear_source(u, T, dx, L=3):
    """Nonconservative terms according to stechmann Majda"""
    from scipy.ndimage import correlate1d

    # compute centered differences
    du = correlate1d(u, [-1/2/dx, 0, 1/2/dx], axis=1)
    dT = correlate1d(T, [-1/2/dx, 0, 1/2/dx], axis=1)

    ft = np.empty_like(T)
    fu = np.empty_like(u)

    fu[0] = 0.0
    fu[1] = 3/2/sqrt(2) * u[1] * du[2]
    fu[2] = 0.0

    ft[0] = 0.0
    ft[1] = -(2 * du[1] * T[2] + T[1] * du[2]/2)/sqrt(2)
    ft[2] = -1/2/sqrt(2) * (u[1] * dT[1]-T[1]*du[1])

    return dict(u=fu, t=ft)


class LinTadmor(Tadmor1DBase):
    def fx(self, out, uc):
        f(uc, fq=out)

class NLTadmor(Tadmor1DBase):

    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def fx(self, out, uc):
        f_nonlinear(uc, fq=out, **self.kwargs)

class MulticloudModel(object):
    """Base MulticloudModel solver class"""
    diags = {}

    def __init__(self, tad=LinTadmor()):
        "docstring"
        self.diags['rms'] = 0.0
        self.reset_avg()
        self.tad = tad

    def update_avg(self, soln):

        if self.avg is None:
            self.avg, self.navg =  soln._data, 1
        else:
            self.avg, self.navg =  (self.avg * self.navg + soln._data)/(self.navg+1), self.navg + 1

    def reset_avg(self):
        self.avg = None
        self.navg = 0.0

    def onestep(self, soln, time, dt, dx, nonlinear=0.0, f=f):
        """Perform a single time step of the multicloud model"""
        from functools import partial

        if not self.validate_soln(soln[:,2:-2]):
            raise ValueError("NAN in solution array")

        # periodic boundaries
        soln.comm()

        # hyperbolic terms
        uc = soln.q[:,:,None]
        uc[:] = self.tad.central_scheme(uc, dx, dt)

        # multicloud model step
        soln['hc'], soln['hd'], soln['lmd'] = mc.multicloud_wrapper(
             soln['fc'], soln['fd'],
             soln['fs'], soln['u'][1],
             soln['u'][2], soln['t'][1],
             soln['t'][2], soln['teb'],
             soln['q'], soln['hs'], dt, dx, time,
             soln['tebst'])

        self.diags['rms'] = np.mean(soln['u']**2)
        self.update_avg(soln)

        return soln

    def init_mc(self, n=1000, dx=40, asst=5, lsst=1.0, **kwargs):

        soln = Soln(n, **kwargs)

        fceq, fdeq, fseq = mc.multicloud_eq()

        logger.info("ASST=`{}` K,n=`{}`, dx=`{}`".format(asst, n, dx))

        asst = asst/15.0
        dx /= 1500

        soln['fc'] = fceq
        soln['fd'] = fdeq
        soln['fs'] = fseq

        # initialize temperature field with small random perturbation
        soln['t'][1] = np.random.randn(n) * .01

        # initialize teb
        x = np.arange(n) * dx
        domain_size = n * dx
        lsst = domain_size * lsst

        tebst = asst * np.cos(2 * np.pi * (x - domain_size / 2) / lsst) * (
            x > domain_size / 2 - lsst / 2) * (x < domain_size / 2 + lsst / 2)
        tebst[x >= domain_size / 2 + lsst / 2] = -asst
        tebst[x <= domain_size / 2 - lsst / 2] = -asst

        soln['tebst'] = tebst

        return soln, dx

    def record_array_soln(self, soln, t):
        return soln.record_array_soln(t)

    def validate_soln(self, soln):
        return not np.any(np.isnan(soln))

class MulticloudModelDissipation(MulticloudModel):
    """Multicloud Model with dissipation"""

    def __init__(self, *args, dissipation=.1, **kwargs):
        "docstring"
        super(MulticloudModelDissipation, self).__init__(*args, **kwargs)

        logger.info("Dissipation=`{}`".format(dissipation))
        self.diags['fdis'] = 0.0
        self.dissipation = dissipation

    def onestep(self, soln, time, dt, dx, *args, **kwargs):
        from functools import partial
        soln = super(MulticloudModelDissipation, self)\
               .onestep(soln, time, dt, dx, *args, **kwargs)

        uold = soln['u'].copy()
        soln['u'] = np.exp(-dt * self.dissipation) * soln['u']
        fdis = (soln['u'] - uold) / dt

        self.diags['fdis'] = np.mean(soln['u'] * fdis)


        return soln



class MulticloudModelNonlinear(MulticloudModelDissipation):
    """Multicloud with Nonlinear advection"""
    def __init__(self, *args, q_tld=.9, **kwargs):
        "docstring"

        tad = NLTadmor(q_tld=q_tld)

        super(MulticloudModelNonlinear, self).__init__(*args, tad=tad, **kwargs)

        # arrays for AB3
        self._nonlinear_ab3 = {key: deque([0,0,0], maxlen=3) for key in ['u', 't']}


        from functools import partial
        logger.info("Q_TLD=`{}`".format(q_tld))
        self._f = partial(f_nonlinear, q_tld=q_tld)

    def _nonlinear_source_update(self, soln, dx, dt):
        soln.comm()
        du = nonlinear_source(soln['u'], soln['t'], dx, L=3)

        for key in du:
            dd = self._nonlinear_ab3[key]
            dd.appendleft(du[key])

            soln[key] = soln[key] + dt * (dd[0] * 23/ 12 - 4/3 * dd[1] + 5/12 *dd[2])

        return soln

    def onestep(self, soln, time, dt, dx, *args, **kwargs):
        soln = super(MulticloudModelNonlinear, self)\
               .onestep(soln, time, dt, dx, f=self._f, *args, **kwargs)

        # compute source terms from nolinear advection
        soln = self._nonlinear_source_update(soln, dx, dt)

        return soln

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


def main(run_duration=100, dt_out=1.0, solver=None, restart_file=None,
         cfl=.1, init_kwargs={}):
    """Runs multicloud model

    TODO This file is too complicated needs to be refactored, and the IO needs
    to be rethought
    """
    t_start = 0.0


    if solver is None:
        solver = MulticloudModel()

    logger.info("solver object = " + repr(solver))
    logger.info("Starting run with duration={0}".format(run_duration))

    soln, dx = solver.init_mc(**init_kwargs)

    if restart_file is not None:
        soln, t_start, dx = load_restart_file(restart_file)
        logger.info("Loading restart file at t={0}".format(t_start))
    elif os.path.exists('ic.npz'):
        soln = solver.init_mc_from_file("ic.npz")

    dt = dx * cfl
    t_end = t_start + run_duration

    t_out = t_start + dt_out

    grid = {
        'x': dx*np.arange(len(soln['teb'])-4),
        'm': np.arange(soln['u'].shape[0])
    }

    # allocate output buffer
    writer = NetCDF4Writer(grid)

    diags = []

    # include initial data
    i_out = 1

    diagfile = open("diags.pkl", "ab")


    for t, soln in steps(solver.onestep, soln, dt, (t_start, t_end), dx):

        solver.diags['time'] = t
        pickle.dump(solver.diags, diagfile)

        # get averages
        if solver.avg is None:
            avg, navg = soln._data, 1
        else:
            avg, navg = (avg * navg + solver.avg)/(navg + 1), navg + 1
        solver.reset_avg()

        if t > t_out:
            writer.collect(t, soln.vars)
            t_out, i_out  = t_out + dt_out, i_out + 1


    writer.collect(t, soln.vars)
    writer.finish()
    save_restart_file("restart_" + str(uuid.uuid1()) + ".pkl", soln, t, dx)
    diagfile.close()

    avg = Soln.from_arr(avg).record_array_soln(t)
    np.savez("averages.npz", avg=avg)



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
