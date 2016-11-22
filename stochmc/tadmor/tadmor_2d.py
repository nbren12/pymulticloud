""" Python implementation of the Tadmor centered scheme in 2d


Routines
--------
central_scheme - 2d implementation of tadmor centered scheme

TODO
----
Replace all periodic_bc calls with `comm', so that this code can be run in parallel
"""
from functools import partial
import numpy as np
from scipy.ndimage import correlate1d
from numba import jit

from .tadmor_common import periodic_bc, minmod


def _roll2d(u):
    return np.roll(np.roll(u, -1, axis=1), -1, axis=2)


@jit
def _slopesy(uc, uy, tht=1.5):
    neq, n1, n2 = uc.shape

    for j in range(neq):
        for i in range(n1):
            for k in range(1, n2 - 1):
                left = tht * (uc[j, i, k + 1] - uc[j, i, k])
                cent = (uc[j, i, k + 1] - uc[j, i, k - 1]) / 2
                right = tht * (uc[j, i, k] - uc[j, i, k - 1])
                uy[j, i, k] = minmod(left, cent, right)


def _slopes(uc, axis=-1, **kwargs):
    uc = np.rollaxis(uc, axis, start=uc.ndim)
    uy = np.empty_like(uc)
    _slopesy(uc, uy, **kwargs)
    return np.rollaxis(uy, -1, start=axis)


def _stagger_avg(uc):
    ux = np.empty_like(uc)
    uy = np.empty_like(uc)

    ux = _slopes(uc, axis=1)
    uy = _slopes(uc, axis=2)

    ox = correlate1d(uc, [.25, .25], axis=1)
    ox += correlate1d(ux, [.125, -.125], axis=1)

    oy = correlate1d(uc, [.25, .25], axis=2)
    oy += correlate1d(uy, [.125, -.125], axis=2)

    return correlate1d(ox, [.5, .5], axis=2) + \
        correlate1d(oy, [.5, .5], axis=1)


def _corrector_step(fx, fy, lmd_x, lmd_y):

    ox = correlate1d(fx, [lmd_x, -lmd_x], axis=1)
    oy = correlate1d(fy, [lmd_y, -lmd_y], axis=2)


    return correlate1d(ox, [.5, .5], axis=2) + \
        correlate1d(oy, [.5, .5], axis=1)


def _roll2d(u):
    return np.roll(np.roll(u, -1, axis=1), -1, axis=2)


class Tadmor2D(object):

    comm = partial(periodic_bc, axes=(1, 2))

    def _single_step(self, fx, fy, uc, dx, dy, dt):
        ux = np.zeros_like(uc)
        uy = np.zeros_like(uc)
        uc = uc.copy()
        lmd_x = dt / dx
        lmd_y = dt / dy

        self.comm(uc)
        ustag = _stagger_avg(uc)

        # predictor: mid-time-step pointewise values at cell-center
        # Eq. (1.1) in Jiand and Tadmor
        ux = _slopes(fx(uc), axis=1)
        uy = _slopes(fy(uc), axis=2)
        uc -= lmd_x / 2 *ux   + lmd_y/2 * uy

        # corrector
        # Eq (1.2) in Jiang and Tadmor
        self.comm(uc)
        ustag += _corrector_step(fx(uc), fy(uc), lmd_x, lmd_y)

        return ustag

    def central_scheme(self, fx, fy, uc, dx, dy, dt):
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
        ustag = _roll2d(self._single_step(fx, fy, uc, dx, dy, dt / 2))
        uc = self._single_step(fx, fy, ustag, dx, dy, dt / 2)

        return uc
