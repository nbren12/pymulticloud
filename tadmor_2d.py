"""
Python implementation of the Tadmor centered scheme in 2d
"""
import numpy as np
from numpy import pi
from scipy.ndimage import correlate1d
import numba
from numba import jit

from tadmor_1d import (periodic_bc, minmod)

def roll2d(u):
    return np.roll(np.roll(u, -1, axis=1), -1,axis=2)



@jit
def slopesy(uc, uy, tht=1.5):
    neq, n1, n2 = uc.shape

    for j in range(neq):
        for i in range(n1):
            for k in range(1, n2 -1):
                left = tht * (uc[j, i, k + 1] - uc[j, i, k])
                cent = (uc[j, i, k + 1] - uc[j, i, k - 1]) / 2
                right = tht * (uc[j, i, k] - uc[j, i, k - 1])
                uy[j, i, k] = minmod(left, cent, right)



def slopes(uc, axis=-1, **kwargs):
    uc = np.rollaxis(uc, axis, start=uc.ndim)
    uy = np.empty_like(uc)
    slopesy(uc, uy, **kwargs)
    return np.rollaxis(uy, -1, start=axis)

def stagger_avg(uc):
    ux = np.empty_like(uc)
    uy = np.empty_like(uc)

    ux = slopes(uc, axis=1)
    uy = slopes(uc, axis=2)

    ox = correlate1d(uc, [.25, .25], axis=1)
    ox += correlate1d(ux, [.125, -.125], axis=1)


    oy = correlate1d(uc, [.25, .25], axis=2)
    oy += correlate1d(uy, [.125, -.125], axis=2)

    return correlate1d(ox, [.5, .5], axis=2) + \
        correlate1d(oy, [.5, .5], axis=1)


def corrector_step(fx, fy, lmd):


    ox = correlate1d(fx, [lmd, -lmd], axis=1)
    oy = correlate1d(fy, [lmd, -lmd], axis=2)


    return correlate1d(ox, [.5, .5], axis=2) + \
        correlate1d(oy, [.5, .5], axis=1)

def single_step(fx, fy, uc, dx, dt):
    ux = np.zeros_like(uc)
    uy = np.zeros_like(uc)
    uc = uc.copy()
    lmd = dt / dx

    periodic_bc(uc, axes=(1, 2))
    ustag = stagger_avg(uc)

    # predictor: mid-time-step pointewise values at cell-center
    # Eq. (1.1) in Jiand and Tadmor
    ux = slopes(fx(uc), axis=1)
    uy = slopes(fy(uc), axis=2)
    uc -= lmd / 2 * (ux + uy)

    # corrector
    # Eq (1.2) in Jiang and Tadmor
    periodic_bc(uc, axes=(1, 2))
    ustag += corrector_step(fx(uc), fy(uc), lmd)

    return ustag

def roll2d(u):
    return np.roll(np.roll(u, -1, axis=1), -1,axis=2)

def central_scheme(fx, fy, uc, dx, dt):
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
    ustag = roll2d(single_step(fx, fy, uc, dx, dt/2))
    uc = single_step(fx, fy, ustag, dx, dt/2)

    return uc

