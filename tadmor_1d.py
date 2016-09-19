"""
Python implementation of the Tadmor centered scheme in 1d
"""
import numpy as np
from numpy import pi 
import numba
from numba import jit

@jit
def periodic_bc(u, g=2):
    u[:, -g:] = u[:, g:2 * g]
    u[:, :g] = u[:, -2 * g:-g]

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
def slopes(uc, ux, tht=2.0):

    neq, n = uc.shape

    for j in range(neq):
        for i in range(1, n - 1):
            left = tht * (uc[j, i + 1] - uc[j, i])
            cent = (uc[j, i + 1] - uc[j, i - 1]) / 2
            right = tht * (uc[j, i] - uc[j, i - 1])

            ux[j, i] = minmod(left, cent, right)

    return ux

def single_step(fx, uc, dx, dt):
    ux = np.zeros_like(uc)
    uc = uc.copy()
    lmd = dt / dx
    periodic_bc(uc)
    slopes(uc, ux)

    # average onto staggered grid
    ustag = ((np.roll(uc, 1, axis=-1) + uc) / 2
             + (ux - np.roll(ux, 1, axis=-1)) / 8)

    # predictor: mid-time-step pointewise values at cell-center
    # Eq. (1.1) in Jiand and Tadmor
    fc = fx(uc)
    slopes(fc, ux)
    uc -= lmd / 2 * ux

    # corrector
    # Eq (1.2) in Jiang and Tadmor
    periodic_bc(uc)
    fc = fx(uc)


    ustag -= lmd * (np.roll(fc, 1, axis=-1) - fc)

    return ustag


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
    ustag =  np.roll(single_step(fx, uc, dx, dt/2), -1 , axis=0)
    uc = single_step(fx, ustag, dx, dt/2)


    return uc




def tadmor_error(n):
    uc = np.zeros((1, n+ 4))
 
    L = 1.0
    dx = L / (n+1)

    x = np.arange(-2,n+2) * dx

    u0 = np.sin(2*pi*x)**10
    uc[0,:] = u0

    dt = dx/2

    def fx(u):
        return u



    for it in range(int(1.0/dt)):
        uc = central_scheme(fx, uc, dx, dt)

    return np.sum(np.abs((uc[0,:]-u0)))/n


def tadmor_convergence():
    import matplotlib.pyplot as plt
    nlist = [500, 1000, 2000, 4000]

    err = [tadmor_error(n) for n in nlist]
    plt.loglog(nlist, err)
    p = np.polyfit(np.log(nlist), np.log( err ), 1)
    plt.title('Order of convergence p = %.2f'%p[0])
    plt.show()



def test_tadmor_1d(n=2000):
    import matplotlib.pyplot as plt
    n = 2000
    uc = np.zeros((1, n+ 4))
 
    L = 1.0
    dx = L / (n+1)

    x = np.arange(-2,n+2) * dx

    uc[0,:] = np.exp(-((x-.5)/.10)**2)

    dt = dx/2


    def fx(u):
        return u


    plt.plot(x, uc[0,:], label='exact')

    for it in range(int(1.0/dt)):
        uc = central_scheme(fx, uc, dx, dt)


    plt.plot(x, uc[0,:], label='tadmor', c='k')

def test_upwind_1d(n=2000):
    import matplotlib.pyplot as plt
    uc = np.zeros((1, n+ 4))
 
    L = 1.0
    dx = L / (n+1)

    x = np.arange(-2,n+2) * dx

    uc[0,:] = np.exp(-((x-.5)/.10)**2)

    dt = dx/2


    def fx(u):
        return u



    for it in range(int(1.0/dt)):
        uc += -dt/dx * (uc - np.roll(uc, -1, axis=-1))


    plt.plot(x, uc[0,:], label='upwind')

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    n = 1000
    test_tadmor_1d(n)
    test_upwind_1d(n)
    plt.legend()
    plt.show()
