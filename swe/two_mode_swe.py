import logging
import numpy as np
from numpy import sqrt
from scipy.ndimage import correlate1d
from .tadmor_1d import periodic_bc, central_scheme

logging.basicConfig(level=logging.INFO)


# include barotropic mode
L = 3


def f(q):
    """Conservative flux according to stechmann and majda"""
    u = q[:2*L:2,...]
    T = q[1:2*L:2,...]

    fq = np.empty_like(q)

    fq[0,...] = 0.0
    fq[1,...] = 0.0
    fq[2,...] = - T[1] + 3/sqrt(2) * u[1] * u[2]
    fq[3,...] = - u[1] + sqrt(2) * u[1] * T[2] - u[2] * T[1] / sqrt(2)
    fq[4,...] = -T[2]
    fq[5,...] = -u[2] / 4


    return fq

def f_nc(q, dx):
    """Nonconservative terms according to stechmann Majda"""
    periodic_bc(q)
    u = q[:2*L:2,...]
    T = q[1:2*L:2,...]

    dq = correlate1d(q, [-1/2/dx, 0, 1/2/dx], axis=1)
    du = dq[:2*L:2,...]
    dT = dq[1:2*L:2,...]

    fq = np.empty_like(q)

    fq[0] = 0.0
    fq[1] = 0.0
    fq[2] = 3/2/sqrt(2) * u[1] * du[2]
    fq[3] = -(2 * du[1] * T[2] + T[1] * du[2]/2)/sqrt(2)
    fq[4] = 0.0
    fq[5] = -1/2/sqrt(2) * (u[1] * dT[1]-T[1]*du[1])

    return fq


def onestep(q, dx, dt):
    q = central_scheme(f, q, dx, dt)
    q += dt * f_nc(q, dx)

    return q

def steps(q, dx, dt, t_end):
    t = 0

    yield t, q

    while (t < t_end - 1e-10):
        dt = min(dt, t_end-t)
        q = onestep(q, dx, dt)
        t+=dt

        logging.info("t = {t:.2f}".format(t=t))

        yield t, q

labs = ['u0','t0', 'u1', 't1', 'u2', 't2']

def main():
    import matplotlib.pyplot as plt
    km = 1000
    Le = 1500 * km
    nx = 500
    dx = 50 * km / Le

    dt = dx * .2
    q = np.zeros((2*L, nx+4))
    x = np.arange(-2, nx+2) * dx
    domain_size = x[-2]
    t_end =  1 * domain_size

    q[3,:] = np.sin(2*np.pi*x/domain_size) * (x < domain_size/2)

    t, q = zip(*list(steps(q, dx, dt, t_end)))

    def plot_q(t, q):
        fig, axs= plt.subplots(2,2)
        for i, ax in enumerate(axs.flat):
            ax.plot(x, q[i + 2] )
            ax.set_title(label=labs[i+2])

        plt.suptitle("t={t}".format(t=t))

    plot_q(t[0], q[0])
    plot_q(t[-1], q[-1])

    plt.figure()
    plt.contourf(np.vstack(qq[3] for qq in q))
    plt.show()

if __name__ == '__main__':
    main()
    pass
