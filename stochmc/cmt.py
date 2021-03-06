"""Implementation of multicloud model with cmt in python

This code uses the forcings and gillespie algorithm from the fortran code

Functions
---------
stochastic_cmt_diagnostic_run(datadir)
stochastic_integrate_array(scmt, rates, a, b):
transition_rates_array(u, qd, qc, lmd)
rhs(t, u, scmt, qd, u_relax):
"""
import sys
import logging
import numpy as np
import pandas as pd
from numpy import log, exp, sqrt, pi, cos
from scipy.interpolate import interp1d
from numba import jit, guvectorize
from random import uniform

logger = logging.getLogger(__file__)

hour = 3600
day = 86400
km = 1000

T = 2.9575e+04
L = 1.4827e+06
c = L / T
alpha_bar = 15.34

tscale = T
qscale = alpha_bar / T

# Useful global variables
hsbar = 0
lmd = 0
dulow = 0

nu = 3

transform_mat = {}
z = {}
nzgrid = 32

z = np.arange(0, nzgrid) * np.pi / nzgrid
m = np.arange(0, nu)

transform_mat[nzgrid] = sqrt(2) * np.cos(z[:, None] * m[None, :])
transform_mat[nzgrid][:, 0] = 1.0


def cos_series_val(z, u):
    z = np.array(z, dtype=np.float64)
    out = np.zeros_like(z)
    for i in range(u.shape[0]):
        out += sqrt(2) * u[i] * cos((i + 1) * z)

    return out


def calculate_dulow(u, a, b, z0=0, **kwargs):
    from scipy.optimize import minimize_scalar

    def f(z):
        return -(cos_series_val(z, u) - cos_series_val(z0, u))**2

    opt = minimize_scalar(f, bounds=(a, b), method='Bounded', **kwargs)
    return opt.x, cos_series_val(opt.x, u) - cos_series_val(z0, u)


@jit
def calc_du(u):

    uz = transform_mat[32].dot(u)
    dulow = np.zeros(u.shape[1])
    dumid = np.zeros(u.shape[1])

    for i in range(uz.shape[1]):
        uzz = uz[:, i]

        iopt = np.abs(uzz[2:14] - uzz[0]).argmax() + 2

        uzst = uzz[iopt]
        dulow[i] = uzst - uzz[0]

        imid = np.abs(uzz[16:] - uzst).argmax() + 16
        dumid[i] = uzz[imid] - uzst

    return dulow, dumid


def test_calculate_dulow():
    import matplotlib.pyplot as plt
    a = 1 * pi / 16
    b = 7 * pi / 16

    u = np.random.rand(3)

    zopt, dulow = calculate_dulow(u, a, b)

    z = np.linspace(a, b, 100)

    plt.plot(cos_series_val(z, u), z)
    plt.plot(cos_series_val(zopt, u), zopt, 'ro')
    plt.show()


@jit(nopython=True)
def heaviside(x):
    return 0.5 * (np.sign(x) + 1)


class rates_decorator(object):
    def __init__(self, m=3, **kwargs):
        "docstring"

        self.jit_kwargs = kwargs
        self.m = m

    def __call__(self, transition_rates):
        """Decorator for jit based transition rates functions
        """
        trates = jit(transition_rates, **self.jit_kwargs)

        @jit
        def func(dulow, qd, qc, lmd, m=self.m):
            n = qd.shape[0]
            rates = np.zeros((m, m, n))

            for i in range(n):
                trates(dulow[i], qd[i], qc[i], lmd[i], rates[:, :, i])

            return rates

        return func


@rates_decorator(nopython=True)
def transition_rates(dulow, qd, qc, lmd, T):
    """Transition rates for stochastic CMT process"""
    taur = 8 * hour
    beta_q = 1 / (10 / day)
    beta_u = 1 / (10.0)
    qcref = 20 / day
    qdref = 20 / day
    duref = 10
    dumin = 5

    dulow = np.abs(dulow)

    T[0, 1] = heaviside(qd) * exp(-lmd + beta_q * qd)
    T[1, 2] = heaviside(dulow - dumin) * exp(beta_u * dulow + beta_q * qc)
    T[1, 0] = exp(lmd + beta_q * (qdref - qd))
    T[2, 0] = T[1, 0]
    T[2, 1] = exp(beta_u * (duref - dulow) + beta_q * (qcref - qc))

    T /= taur


@jit
def softmax(x):
    return (np.abs(x) + x) / 2


@rates_decorator(nopython=True)
def tr1(dulow, qd, qc, lmd, T):
    """Transition rates for stochastic CMT process

    This function is empirically hand-tuned. Written on 2016-10-17.
    """
    q = qd + .4 * qc

    a = 20 * softmax(lmd + .5)
    b = np.tanh(softmax(q - 1))
    bn = np.tanh(softmax(1 - q))

    c = softmax(abs(dulow) - .3) * 10
    cn = softmax(.1 - abs(dulow)) * 10

    T[0, 1] = a + b
    T[1, 2] = c
    T[1, 0] = bn
    T[2, 0] = bn
    T[2, 1] = cn

    T /= day


@rates_decorator(nopython=True)
def tr2(dulow, qd, qc, lmd, T):
    """Transition rates for stochastic CMT process

    This function is empirically hand-tuned. Written on 2016-10-18.
    """
    q = qd + .4 * qc

    a = 20 * softmax(lmd + .6)
    b = np.tanh(softmax(q - .5))
    bn = np.tanh(softmax(.5 - q))

    c = softmax(abs(dulow) - .15) * 15
    cn = softmax(.05 - abs(dulow)) * 15

    T[0, 1] = a + b
    T[1, 2] = c
    T[1, 0] = bn
    T[2, 0] = bn
    T[2, 1] = cn

    T /= day


@jit(nopython=True)
def update_cmt(u, scmt, qd, dulow, dumid, dt):

    d1 = 1 / (3 * day)
    d2 = 1 / (3 * day)
    tauf = 1.25 * day
    qdref = 10 / day

    for j in range(u.shape[0]):
        for i in range(u.shape[1]):

            if scmt[i] == 0:
                u[j, i] = exp(-d1 * dt) * u[j, i]
            elif scmt[i] == 1:
                u[j, i] = exp(-d2 * dt) * u[j, i]

            elif scmt[i] == 2:
                if dumid[i] * dulow[i] < 0:
                    kappa = -np.tanh((qd[i] / qdref)**2) * dumid[i] / tauf
                else:
                    kappa = 0

                u[1, i] = u[1, i] + kappa * dt

    return u


def stochastic_integrate_array(scmt, rates, a, b):
    """Preform Gillespie algorithm for an array of data

    Parameters
    ----------
    scmt: (n,)
        array of integers representing the discrete stochastic state
    rates: (p, p, n)
        array of stochastic transition rates. p is the number of allowed states
    a: float
        starting time
    b: float
        ending time

    Returns
    -------
    scmt: (n,)
       updated stochastic state array
    """

    n = scmt.shape[0]
    time = np.ones(n) * a
    running = np.arange(scmt.shape[0])

    crates = np.cumsum(rates, axis=1)

    while len(running) > 0:
        logger.debug("running has length {0}".format(len(running)))

        lam = crates[scmt[running], :, np.arange(len(running))]

        mask = lam[:, -1] != 0

        running, lam = running[mask], lam[mask, :]

        U1 = np.random.rand(len(running))
        tau = -log(U1) / lam[:, -1]
        time[running] = time[running] + tau

        mask = time[running] < b
        running = running[mask]

        if len(running) > 0:
            lam = lam[mask, :]

            U2 = np.random.rand(len(running))

            for i, idx in enumerate(running):
                scmt[idx] = np.searchsorted(lam[i, :], U2[i] * lam[i, -1])

    return scmt

    # This version of the code goes grid point by grid point
    # for i in range(scmt.shape[0]):
    #     t = a
    #     transition_rates(aulow[i], qd[i], qc[i], lmd[i], rates)
    #     cum_rates = np.cumsum(rates, axis=1)
    #     while True:
    #         cs = cum_rates[scmt[i],:]
    #         U = uniform(0.0, 1.0)
    #         tau = -log(U) / cs[-1]

    #         if t + tau < b:
    #             t += tau
    #             U = uniform(0.0,1.0)
    #             action_index = np.searchsorted(cs, U * cs[-1])
    #             scmt[i] = action_index
    #         else:
    #             break

    # return scmt


def stochastic_integrate(scmt, dulow, qd, qc, lmd, a, b):
    """Stochastic integration using gillespie algorithm"""

    # stochastic integration
    # uses gillespie's algorithm
    t = a

    # do while loop
    while True:
        rates = transition_rates(dulow, qd, qc, lmd)
        rates = np.cumsum(rates[scmt, :])
        U = uniform(0.0, 1.0)
        tau = -log(U) / rates[-1]

        if t + tau < b:
            t += tau
            U = uniform(0.0, 1.0)
            action_index = np.searchsorted(rates, U * rates[-1])
            scmt = action_index
        else:
            break

    return scmt


def interpolant_hash(tout, qd, dt_in):
    """ Create a hash table of the needed Qd values
    """
    tout_iter = tout.flat
    t = next(tout_iter)
    cache_times = [t]
    for i, next_time in enumerate(tout_iter):
        while t < next_time - 1e-10:
            dt = min(next_time - t, dt_in)
            t += dt
            cache_times.append(t)

    qd_interp = qd(np.array(cache_times))
    qd_cache = dict(zip(cache_times, qd_interp))

    return qd_cache


def run_cmt_model(u, scmt, tout, *args, f=None, dt_in=600):

    t0 = tout[0]

    if scmt.shape != u(t0).shape[1:]:
        raise ValueError("SCMT must have the same shape as u")

    scmt = scmt.astype(np.int32)


    output_scmt = np.zeros((tout.shape[0], u(t0).shape[1]))

    # Precalculate qd at necessary times
    caches = [interpolant_hash(tout, arg, dt_in) for arg in args]
    u_cache = interpolant_hash(tout, u, dt_in)

    if f is None:
        f = tr1

    tout_iter = tout.flat
    t = next(tout_iter)
    for i, next_time in enumerate(tout_iter):
        while t < next_time - 1e-10:
            dt = min(next_time - t, dt_in)

            arg_t = [cache[t] for cache in caches]
            ut = u_cache[t]

            dulow, duhi = calc_du(ut)

            # stochastic integration
            rates = f(dulow, *arg_t)
            scmt = stochastic_integrate_array(scmt, rates, t, t + dt)
            t += dt

        # store output
        logger.info("Time %.2f days" % (next_time / day))
        output_scmt[i] = scmt

    return output_scmt


def stochastic_cmt_diagnostic_run(data, f=tr1):
    """Run stochastic cmt scheme in diagnostic mode on previous output from the
    multicloud model.

    Parameters
    ----------
    datadir: str
        directory containing output from the python version of the multicloud model
    """
    from .read import read_data
    import matplotlib.pyplot as plt
    import xarray as xr

    # read qd
    logger.info("starting column cmt run")

    t_data = data.time.values * T
    qd_data = data.hd.values
    qc_data = data.hc.values
    lmd_data = data.lmd.values

    u_relax_data = data['u']

    qc = interp1d(t_data, qc_data, axis=0)
    qd = interp1d(t_data, qd_data, axis=0)
    u_relax = interp1d(t_data, u_relax_data, axis=0)
    lmd = interp1d(t_data, lmd_data, axis=0)

    # initialization
    n = qd_data.shape[1]
    scmt = np.zeros(n, dtype=np.int32)

    # tout
    tout = t_data
    output_cmt = run_cmt_model(u_relax, scmt, tout, qd, qc, lmd, dt_in=600, f=f)

    return tout, output_cmt


class CmtSolver(object):
    def __init__(self, transition_rates=tr2, mcsolver=None):
        "docstring"

        logger.info("Using CMT transition rate function: " + repr(
            transition_rates.__code__))

        if mcsolver is None:
            from .swe.multicloud import MulticloudModel
            mcsolver = MulticloudModel()
        logger.info("Using MulticloudSolver "+ repr(mcsolver))

        self._multicloud_model = mcsolver
        self._du = [0, 0, 0]
        self.diags['kcmt'] = np.zeros((3, ))
        self.diags['ncmt'] = np.zeros((3, ), dtype=np.int32)

        self.transition_rates = transition_rates

    def init_mc(self, *args, **kwargs):

        soln, dx = self._multicloud_model.init_mc(*args,
                                                  extra_vars=['scmt'],
                                                  **kwargs)

        return soln, dx

    def onestep(self, soln, time, dt, *args, **kwargs):

        soln = self._multicloud_model.onestep(soln, time, dt, *args, **kwargs)
        soln = self._cmt_step(soln, time, dt)

        return soln

    def _cmt_step(self, soln, time, dt):
        """Step of cmt model"""

        u = soln['u']
        hd = soln['hd']
        hc = soln['hc']
        lmd = soln['lmd']
        scmt = soln['scmt'].astype(np.int32)

        dulow, dumid = calc_du(u)
        rates = self.transition_rates(dulow, hd, hc, lmd)
        scmt = stochastic_integrate_array(scmt, rates, T * time,
                                          T * (time + dt))

        uold = u.copy()

        u = update_cmt(u, scmt, hd, dulow, dumid, dt * T)

        # diagnostic

        fcmt = (u - uold) / dt


        # compute cmt diagnostics
        for i in range(3):
            mask = scmt == i
            nmask = mask.sum()
            self.diags['ncmt'][i] = nmask
            if nmask > 0:
                self.diags['kcmt'][i] = np.mean(fcmt[:, mask] * u[:, mask])
            else:
                self.diags['kcmt'][i] = 0.0


        soln['u'] = u
        soln['scmt'] = scmt.astype(np.float64)

        return soln

    def __getattr__(self, name):
        return self._multicloud_model.__getattribute__(name)


def main():
    t, scmt = stochastic_cmt_diagnostic_run("data")
    np.savez("scmt.npz", t=t, scmt=scmt)


if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO)
    main()
