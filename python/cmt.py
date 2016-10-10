"""Implementation of multicloud model with cmt in python

This code uses the forcings and gillespie algorithm from the fortran code
"""
import sys
import logging
import numpy as np
import pandas as pd
from numpy import log, exp, sqrt, pi, cos
from scipy.interpolate import interp1d
from numba import jit
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
qscale = alpha_bar

# Useful global variables
hsbar = 0
lmd = 0
dulow = 0

nu = 3

transform_mat = {}; z={};
nzgrid = 32

z = np.arange(0,nzgrid) * np.pi/nzgrid
m = np.arange(0, nu)

transform_mat[nzgrid] = sqrt(2) * np.cos(z[:,None] * m[None,:])
transform_mat[nzgrid][:,0] = 1.0



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
        uzz = uz[:,i]

        iopt = np.abs(uzz[2:14] - uzz[0]).argmax()

        uzst = uzz[iopt]
        dulow[i] = uzst-uzz[0]


        imid = np.abs(uzz[16:]-uzst).argmax()
        dumid[i] = uzz[imid]-uzz[iopt]

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

@jit(nopython=True)
def transition_rates(dulow, qd, qc, lmd, T):
    """Transition rates for stochastic CMT process"""
    taur = 8 * hour
    beta_lmd = 1
    beta_q = 1 / (10 / day)
    beta_u = 1 / (10.0)
    qcref = 10 / day
    qdref = 10 / day
    duref = 20
    dumin = 5


    dulow = abs(dulow)

    T[0, 1] = heaviside(qd) * exp(beta_lmd * (1 - lmd) + beta_q * qd)
    T[1, 2] = heaviside(dulow - dumin) * exp(beta_u * dulow + beta_q * qc)

    T[1, 0] = exp(beta_lmd * lmd + beta_q * (qdref - qd))
    T[2, 0] = T[1, 0]
    T[2, 1] = exp(beta_u * (duref - dulow) + beta_q * (qcref - qc))

    T /= taur

@jit
def transition_rates_array(u, qd, qc, lmd):
    n = qd.shape[0]
    rates = np.zeros((3,3, n))
    dulow, dumid = calc_du(u)

    for i in range(n):
        transition_rates(dulow[i], qd[i], qc[i], lmd[i], rates[:,:,i])

    return rates

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
        dulow, dumid = calc_du(u)

        du = np.zeros_like(u)
        if dumid * dulow < 0:
            kappa = -(qd/ qdref)**2 * dumid / tauf

            du[0] = kappa
            du[1] = 0.0
            du[2] = -kappa
    else:
        raise ValueError("scmt must be either 0, 1, or 2")

    return du

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
    """

    n = scmt.shape[0]
    time = np.ones(n) * a
    running = np.arange(scmt.shape[0])

    crates = np.cumsum(rates, axis=1)

    while len(running) > 0:
        logger.debug("running has length {0}".format(len(running)))

        lam = crates[scmt[running],
                    :,
                    np.arange(len(running))]

        mask = lam[:,-1] != 0

        running = running[mask]
        lam = lam[mask, :]

        U1 = np.random.rand(len(running))
        tau = -log(U1)/lam[:,-1]
        time[running] = time[running] + tau

        mask = time[running] < b
        running = running[mask]

        if len(running) > 0:
            lam = lam[mask,:]

            U2 = np.random.rand(len(running))

            for i, idx in enumerate(running):
                scmt[idx] = np.searchsorted(lam[i,:], U2[i]*lam[i,-1])


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
        U = uniform(0.0,1.0)
        tau = -log(U) / rates[-1]

        if t + tau < b:
            t += tau
            U = uniform(0.0,1.0)
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


    output_scmt = np.zeros((tout.shape[0], u(0).shape[1]))

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
            rates = transition_rates_array(ut, qdt, qct, lmdt)
            scmt = stochastic_integrate_array(scmt, rates, t, t + dt)
            t += dt

        # store output
        logger.info("Time %.2f days"%(next_time/day))
        output_scmt[i] = scmt

    return output_scmt


def main_stochastic_range(datadir):
    from .read import read_data
    import matplotlib.pyplot as plt


    # read qd
    logger.info("starting column cmt run")
    data = read_data(datadir)

    t_data = data['time'][:] * T
    qd_data = data['hd'] * alpha_bar / T
    qc_data = data['hc'] * alpha_bar / T
    lmd_data = .4 * np.ones_like(qd_data)
    u_relax_data = data['u'] * c


    qc = interp1d(t_data, qc_data, axis=0)
    qd = interp1d(t_data, qd_data, axis=0)
    u_relax = interp1d(t_data, u_relax_data, axis=0)
    lmd =  interp1d(t_data, lmd_data, axis=0)




    # initialization
    n = qd_data.shape[1]
    scmt = np.zeros(n, dtype=np.int32)

    # tout
    tout = np.mgrid[t_data.min():t_data.max():hour]
    output_cmt = run_cmt_model(u_relax, scmt, tout, qd, qc, lmd, dt_in=600)

    return tout, output_cmt



if __name__ == '__main__':
    # test_calculate_dulow()
    logging.basicConfig(level=logging.INFO)
    cProfile.run("""
main_stochastic_range("data/")
        """, "cmt.prof")

#     pdb.run("""
# main_stochastic_range("data/")
#         """)
