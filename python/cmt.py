"""Column model for Majda and Stechmann (2008)


Note
----
Does not need interactive heating from multicloud model.
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
            kappa = -(qd / qdref)**2 * dumid / tauf

            du[0] = kappa
            du[1] = 0.0
            du[2] = -kappa
    else:
        raise ValueError("scmt must be either 0, 1, or 2")

    return du


def test_transition_rates():
    T = transition_rates(10, 10 / day, 10 / day, .4)


def plot_qdtimeseries():
    import matplotlib.pyplot as plt
    qd = QdTimeSeries()

    t = np.mgrid[0:100 * day:60]

    plt.plot(t, qd(t))
    plt.show()


logger.info("Generating QdTimeSeries")


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


def run_column_model(u, scmt, tout, qd, qc, lmd, u_relax=None, dt_in=600):

    # relax to intial profile
    if u_relax is None:
        u_relax = u.copy()

    output = np.zeros((tout.shape[0], len(u)))
    output_scmt = np.zeros((tout.shape[0], ))
    output[0, :] = u

    # Precalculate qd at necessary times
    qd_cache = interpolant_hash(tout, qd, dt_in)
    qc_cache = interpolant_hash(tout, qc, dt_in)
    lmd_cache = interpolant_hash(tout, lmd, dt_in)

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
            fu = rhs(t, u, scmt, qdt, u_relax)
            u += dt * fu

            # store output
        logger.info("Time %.2f days" % (next_time / day))
        output[i, :] = u
        output_scmt[i] = scmt

    return output, output_scmt


def output_to_dataframe(t, output, output_cmt):

    df = pd.DataFrame(
        {'u1': output[:, 0],
         'u2': output[:, 1],
         'u3': output[:, 2],
         'cmt': output_cmt},
        index=t)
    return df


def main_stochastic(heating_filename, output_filename):
    import matplotlib.pyplot as plt

    # read qd
    logger.info("Reading in and interpolating Qd time series")
    d = np.load(heating_filename)
    qd = interp1d(d['t'], d['qd'])

    if 'qc' in d:
        qc = interp1d(d['t'], d['qc'])
    else:
        qc = qd

    if 'lmd' in d:
        lmd = interp1d(d['t'], d['lmd'])
    else:
        lmd = lambda t: 0.0 * t + 0.4

    # initialization
    u = np.zeros((3, ))
    u[0] = 10
    u[1] = -10

    scmt = 0

    # tout
    tout = np.mgrid[d['t'].min():d['t'].max():hour]

    logger.info("Starting column model run")
    output, output_cmt = run_column_model(u,
                                          scmt,
                                          tout,
                                          qd,
                                          qc,
                                          lmd,
                                          dt_in=600)

    output_to_dataframe(tout, output, output_cmt).to_csv(output_filename)


if __name__ == '__main__':
    # test_calculate_dulow()
    main_stochastic(sys.argv[1], sys.argv[2])
