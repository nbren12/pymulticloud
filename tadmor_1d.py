"""
Python implementation of the Tadmor centered scheme in 1d
"""
import numpy as np
from numpy import pi
from scipy.ndimage import correlate1d
import numba
from numba import jit

def periodic_bc(u, g=2, axes=(1,)):
    """periodic bc in arbitrary dimensions

    TODO
    ----
    - refactor this function into another module

    """

    for i in axes:
        idx_in = [slice(None)]*u.ndim
        idx_out = [slice(None)]*u.ndim

        idx_in[i] = slice(g,2*g)
        idx_out[i] = slice(-g, None)

        u[idx_out] = u[idx_in]

        idx_in[i] = slice(-2*g,-g)
        idx_out[i] = slice(0, g)
        u[idx_out] = u[idx_in]
@jit
def minmod(*args):
    mmin = min(args)
    mmax = max(args)

    if mmin > 0.0:
        return mmin
    elif mmax < 0.0:
        return mmax
    else:
        return 0.0

@jit
def _slopes(uc, ux, tht=1.5):

    neq, n = uc.shape

    for j in range(neq):
        for i in range(1, n - 1):
            left = tht * (uc[j, i + 1] - uc[j, i])
            cent = (uc[j, i + 1] - uc[j, i - 1]) / 2
            right = tht * (uc[j, i] - uc[j, i - 1])

            ux[j, i] = minmod(left, cent, right)

    return ux

def _stagger_avg(uc):
    ux = np.empty_like(uc)
    _slopes(uc, ux)
    ustag = (correlate1d(uc,[.5, .5], origin=0, axis=1) +
             correlate1d(ux, [.125, -.125], origin=0, axis=1))

    return ustag

def _cent_avg(ustag):
    """Inverse operation of _stagger_avg"""
    return np.roll(_stagger_avg(ustag), -1, axis=-1)


def _single_step(fx, uc, dx, dt):
    ux = np.zeros_like(uc)
    uc = uc.copy()
    lmd = dt / dx

    periodic_bc(uc)
    ustag = _stagger_avg(uc)


    # predictor: mid-time-step pointewise values at cell-center
    # Eq. (1.1) in Jiand and Tadmor
    fc = fx(uc)
    _slopes(fc, ux)
    uc -= lmd / 2 * ux

    # corrector
    # Eq (1.2) in Jiang and Tadmor
    periodic_bc(uc)
    fc = fx(uc)
    ustag -= lmd * (fc - np.roll(fc, 1, axis=-1))
    return ustag

    # periodic_bc(ustag)
    # uc = _cent_avg(ustag)

    # return uc




def central_scheme(fx, uc, dx, dt):
    """ One timestep of centered scheme


    Parameters
    ----------
    fx : callable
        fx(u) calculates the numeric flux in the x-direction
    uc: (neq, n)
        The state vector on the centered grid
    dx: float
        size of grid cell
    dt: float
        Time step

    Returns
    -------
    out: (neq, n)
       state vector on centered grid
    """
    ustag = np.roll(_single_step(fx, uc, dx, dt/2), -1, axis=-1)
    uc = _single_step(fx, ustag, dx, dt/2)

    # uc = _single_step(fx, uc, dx, dt)

    return uc


def tadmor_error(n):
    uc = np.zeros((1, n+ 4))

    L = 1.0
    dx = L/n

    x = np.arange(-2,n+2) * dx

    u0 = np.sin(2*pi*x)**10
    uc[0,:] = u0

    dt = dx/2

    def fx(u):
        return u

    tend = .87
    t = 0
    while (t < tend - 1e-10):
        dt = min(dt, tend-t)
        uc = central_scheme(fx, uc, dx, dt)
        t+=dt

    uexact = np.sin(2*pi*(x-t))**10

    return np.sum(np.abs((uc[0,:]-uexact)))/n


def test_tadmor_convergence(plot=False):
    """
    Create error convergence plots for 1d advection problem
    """
    nlist = [50, 100, 200, 400, 800, 1600]

    err = [tadmor_error(n) for n in nlist]
    p = np.polyfit(np.log(nlist), np.log( err ), 1)

    if - p[0] < 1.9:
        raise ValueError('Order of convergence is less than 2')
    if plot:
        import matplotlib.pyplot as plt
        plt.loglog(nlist, err)
        plt.title('Order of convergence p = %.2f'%p[0])
        plt.show()



def plot_tadmor_1d(n=2000):
    """
    scalar advection for tadmor scheme
    """
    import matplotlib.pyplot as plt
    uc = np.zeros((1, n+ 4))


    L = 1.0
    dx = L /n

    x = np.arange(-2,n+2) * dx

    uc[0,:] = np.exp(-((x-.5)/.10)**2)

    dt = dx/2


    def fx(u):
        return u


    plt.plot(x, uc[0,:], label='exact')

    tend = 1.8
    t = 0
    while (t < tend - 1e-10):
        dt = min(dt, tend-t)
        uc = central_scheme(fx, uc, dx, dt)
        t+=dt

    plt.plot(x, uc[0,:], label='tadmor', c='k')

def plot_upwind_1d(n=2000):
    """
    scalar advection for upwinding scheme
    """
    import matplotlib.pyplot as plt
    uc = np.zeros((1, n+ 4))

    L = 1.0
    dx = L/n

    x = np.arange(-2,n+2) * dx

    uc[0,:] = np.exp(-((x-.5)/.10)**2)

    dt = dx/2

    tend = 1.8
    t = 0
    while (t < tend - 1e-16):
        dt = min(tend-t, dt)
        uc += dt/dx * (np.roll(uc, 1, axis=-1) - uc)
        t += dt

    plt.plot(x, uc[0,:], label='upwind')

def compare_upwind_tadmor():
    import matplotlib.pyplot as plt

    n = 200
    plot_tadmor_1d(n)
    plot_upwind_1d(n)
    plt.legend()
    plt.show()
