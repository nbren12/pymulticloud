"""Implementation of multicloud model with cmt in python

This code uses the forcings and gillespie algorithm from the fortran code
"""
import sys
import logging
import numpy as np
import pandas as pd
from numpy import log, exp, sqrt, pi, cos
from scipy.interpolate import interp1d

logger = logging.getLogger(__file__)

hour = 3600
day = 86400
km = 1000

T = 2.9575e+04
L = 1.4827e+06
c = L / T
alpha_bar = 15.34

tscale = T
qscale = alpha_bar

# Useful global variables
hsbar = 0
lmd = 0
dulow = 0

nu = 3


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


def calculate_dulow_approx(u, a, b, z0=0, **kwargs):

    def f(z):
        return -(cos_series_val(z, u) - cos_series_val(z0, u))**2

    zgrid = np.linspace(a, b, 8)

    iopt = f(zgrid).argmin()
    zopt = zgrid[iopt]

    return zopt, cos_series_val(zopt, u) - cos_series_val(z0, u)


def test_calculate_dulow():
    import matplotlib.pyplot as plt
    a = 1 * pi / 16
    b = 7 * pi / 16

    u = np.random.rand(3)

    zopt, dulow = calculate_dulow(u, a, b)
    # zopt, dulow = calculate_dulow_approx(u, a, b)

    z = np.linspace(a, b, 100)

    plt.plot(cos_series_val(z, u), z)
    plt.plot(cos_series_val(zopt, u), zopt, 'ro')
    plt.show()

def heaviside(x):
    return 0.5 * (np.sign(x) + 1)


def transition_rates(dulow, qd, qc, lmd):
    """Transition rates for stochastic CMT process"""
    taur = 8 * hour
    beta_lmd = 1
    beta_q = 1 / (10 / day)
    beta_u = 1 / (10.0)
    qcref = 10 / day
    qdref = 10 / day
    duref = 20
    dumin = 5

    T = np.zeros((3, 3))

    dulow = abs(dulow)

    T[0, 1] = heaviside(qd) * exp(beta_lmd * (1 - lmd) + beta_q * qd)
    T[1, 2] = heaviside(dulow - dumin) * exp(beta_u * dulow + beta_q * qc)

    T[1, 0] = exp(beta_lmd * lmd + beta_q * (qdref - qd))
    T[2, 0] = T[1, 0]
    T[2, 1] = exp(beta_u * (duref - dulow) + beta_q * (qcref - qc))

    T /= taur

    return T


def rhs(t, u, scmt, qd, u_relax):

    d1 = 1 / (3 * day)
    d2 = 1 / (3 * day)
    tauf = 1.25 * day
    qdref = 10 / day

    if scmt == 0:
        du = -d1 * (u - u_relax)
    elif scmt == 1:
        du = -d2 * (u - u_relax)
    elif scmt == 2:
        # calculate dumid
        zst, dulow = calculate_dulow_approx(u, pi / 16, 7 * pi / 16, z0=0)
        _, dumid = calculate_dulow_approx(u, pi / 16 * 7, pi / 16 * 13, z0=zst)

        du = np.zeros_like(u)
        if dumid * dulow < 0:
            kappa = -(qd/ qdref)**2 * dumid / tauf

            du[0] = kappa
            du[1] = 0.0
            du[2] = -kappa
    else:
        raise ValueError("scmt must be either 0, 1, or 2")

    return du


def stochastic_integrate(scmt, dulow, a, b, qd, qc, lmd):
    """Stochastic integration using gillespie algorithm"""
    from numpy.random import rand

    # stochastic integration
    # uses gillespie's algorithm
    t = a

    # do while loop
    while True:
        rates = transition_rates(dulow, qd, qc, lmd)
        rates = np.cumsum(rates[scmt, :])
        U = rand()
        tau = -log(U) / rates[-1]

        if t + tau < b:
            t += tau
            U = rand()
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



def run_cmt_model(u, scmt, tout, qd, qc, lmd, dt_in=600):

    output_scmt = np.zeros((tout.shape[0], ))

    # Precalculate qd at necessary times
    qd_cache = interpolant_hash(tout, qd, dt_in)
    qc_cache = interpolant_hash(tout, qc, dt_in)
    lmd_cache = interpolant_hash(tout, lmd, dt_in)
    u_cache = interpolant_hash(tout, u, dt_in)


    tout_iter = tout.flat
    t = next(tout_iter)
    for i, next_time in enumerate(tout_iter):
        while t < next_time - 1e-10:
            dt = min(next_time - t, dt_in)

            qdt = qd_cache[t]
            qct = qc_cache[t]
            lmdt = lmd_cache[t]
            ut  = u_cache[t]

            # stochastic integration
            zst, dulow = calculate_dulow_approx(ut, pi / 16, 7 * pi / 16, z0=0)
            scmt = stochastic_integrate(scmt, dulow, t, t + dt, qdt, qct, lmdt)
            t += dt

        # store output
        logger.info("Time %.2f days"%(next_time/day))
        output_scmt[i] = scmt

    return output_scmt

def run_column_model(u, scmt, tout, qd, qc, lmd, u_relax=None, dt_in=600):

    output = np.zeros((tout.shape[0], len(u)))
    output_scmt = np.zeros((tout.shape[0], ))
    output[0, :] = u

    # Precalculate qd at necessary times
    qd_cache = interpolant_hash(tout, qd, dt_in)
    qc_cache = interpolant_hash(tout, qc, dt_in)
    lmd_cache = interpolant_hash(tout, lmd, dt_in)
    u_relax_cache = interpolant_hash(tout, u_relax, dt_in)


    tout_iter = tout.flat
    t = next(tout_iter)
    for i, next_time in enumerate(tout_iter):
        while t < next_time - 1e-10:
            dt = min(next_time - t, dt_in)


            qdt = qd_cache[t]
            qct = qc_cache[t]
            lmdt = lmd_cache[t]

            # stochastic integration
            zst, dulow = calculate_dulow_approx(u, pi / 16, 7 * pi / 16, z0=0)
            scmt = stochastic_integrate(scmt, dulow, t, t + dt, qdt, qct, lmdt)
            t += dt

            # deterministic integration
            fu = rhs(t, u, scmt, qdt, u_relax_cache[t])
            u += dt * fu

            # store output
        logger.info("Time %.2f days"%(next_time/day))
        output[i, :] = u
        output_scmt[i] = scmt

    return output, output_scmt

def output_to_dataframe(t, output, output_cmt):
    
    df = pd.DataFrame({'u1': output[:,0],
                  'u2': output[:,1],
                  'u3': output[:,2],
                  'cmt': output_cmt}, index=t)
    return df


def main_stochastic(datadir, iloc=500, stochonly=True):
    from .read import read_data
    import matplotlib.pyplot as plt


    # read qd
    logger.info("starting column cmt run for iloc={0}".format(iloc))
    data = read_data(datadir)

    t_data = data['time'][:] * T
    qd_data = data['hd'][:,iloc] * alpha_bar / T
    qc_data = data['hc'][:,iloc] * alpha_bar / T

    u_relax_data = np.zeros((t_data.shape[0], 3))
    u_relax_data[:,:2] = data['u'][:,1:3,iloc] * c


    qc = interp1d(t_data, qc_data)
    qd = interp1d(t_data, qd_data)
    u_relax = interp1d(t_data, u_relax_data, axis=0)

    lmd = lambda t: 0.0 * t + 0.4

    # initialization
    u = np.zeros((3, ))
    u[0] = 0
    u[1] = 0

    scmt = 0

    # tout
    tout = np.mgrid[t_data.min():t_data.max():hour]

    if stochonly:
        output_cmt = run_cmt_model(u_relax, scmt, tout, qd, qc, lmd, dt_in=600)
        df = pd.DataFrame({'cmt': output_cmt}, index=tout)
    else:
        logger.info("Starting column model run")
        output, output_cmt = run_column_model(u, scmt, tout, qd, qc, lmd,
                                            u_relax=u_relax, dt_in=600)

        df = output_to_dataframe(tout, output, output_cmt)

    return df



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    # test_calculate_dulow()
    main_stochastic("data/")
