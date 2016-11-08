"""This module contains solvers for the two mode SWE equations presented in
Stechman and Majda

The barotropic pressure will be solved using the projection method of Chorin.
The poisson solve is done in spectral space. This module assumes that all
arrays have the shape (...,ny, nx). In other words, the final axis of the array
is the x direction. The x-direction is periodic and the y direction has a wall.

"""

from math import pi
import numpy as np
from numba import jit
from .tadmor.tadmor_2d import central_scheme
from scipy.fftpack import dst, idst, dct, idct, fft, ifft, fftfreq



def divergence(u,v, scalex=1.0, scaley=1.0):


    fv = xform_sin(v)
    fu = xform_cos(u)

    k = fftfreq(fv.shape[1], 1/fv.shape[1])[None,:]/scalex
    m = np.arange(fv.shape[0])[:,None]/scaley

    fd = fv * m + fu * 1j * k

    return ixform_cos(fd)

def xform_sin(v):
    fv = np.zeros_like(v, dtype=np.complex128)
    fv[1:-1,:] = idst(fft(v[1:-1,:], axis=1), axis=0, type=1)
    
    return fv

def ixform_sin(v):
    fv = np.zeros_like(v, dtype=np.float64)
    fv[1:-1,:] = dst(np.real(ifft(v[1:-1,:], axis=1)), axis=0, type=1)/2/(fv.shape[0]-1)
    
    return fv

def xform_cos(u):
    return  idct(fft(u, axis=1), type=1, axis=0)


def ixform_cos(u):
    return  dct(ifft(u, axis=1), type=1, axis=0)/2/(u.shape[0]-1)


def pressure_solve_2d(u, v, scalex=1.0, scaley=1.0):
    """This solver assumes a rigid wall at y=+- Ly/2 and periodic boundary
    conditions in the x-direction. Currently assumes that x in [0, 2pi*scalyx] and y
    in [0, pi*scaley].

    Parameters
    ----------
    u: [ny, nx]
        zonal velocity.
    v: [ny, nx]
        meridional velocity. must be zero at y-boundaries.
    scalex: float
    scaley: float

    Returns
    -------
    px: horizontal pressure gradient
    py: horizontal pressure gradient
    p:  pressure
    """


    fv = xform_sin(v)
    fu = xform_cos(u)

    k = fftfreq(fv.shape[1], 1/fv.shape[1])[None,:]/scalex
    m = np.arange(fv.shape[0])[:,None]/scaley

    lapl = -(k**2 + m**2)
    lapl[0,0] = 1.0

    fp = (1j * k *fu + m * fv)/lapl
    fpx = (-k**2 *fu + 1j * k * m * fv)/lapl
    fpy = (-1j * k * m * fu - m**2 * fv )/lapl

    p  = ixform_cos(fp)
    px = ixform_cos(fpx)
    py = ixform_sin(fpy)

    return px, py, p


def Grid(object):
    """Object for transparent handling of pressure solving
    """

    def __init__(self, scales=None, size=None):
        self.scales = scales
        self.size = size

    @property
    def grid(self):
        """
        y, x tuple of meshgrid output
        """
        ny, nx = self.size

        x  = np.arange(nx)/nx * 2 * pi * self.scales[1]
        y  = np.linspace(0, pi, ny) * self.scales[0]


        return np.meshgrid(y, x, indexing='ij')

    def pressure_solve_2d(self, u, v):
        return pressure_solve_2d(u, v, scalex=self.scales[1],
                                 scaley=self.scales[0])

    def divergence(self, u, v):
        return divergence(u, v, scalex=self.scales[1],
                                 scaley=self.scales[0])


def test_pressure_solve_2d():
    x  = np.arange(102)/102 *2* pi
    y = np.linspace(0,pi,51)


    y, x= np.meshgrid(y, x, indexing='ij')


    v = np.sin(4*y) * np.cos(5*x)
    u = 0.0*v # np.cos(3*y) * np.cos(2*x)


    pxex = 4 * np.cos(4*y) *5* np.sin(5*x)/(16+25)
    pyex = 16 * np.sin(4*y) * np.cos(5*x)/(16+25)

    px, py,_ = pressure_solve_2d(u, v, 1.0, 1.0)

    ut, vt = u-px, v-py
    d1 = divergence(ut, vt)

    np.testing.assert_almost_equal(pxex, px)
    np.testing.assert_almost_equal(pyex, py)
    np.testing.assert_almost_equal(d1, 0.0)


    # test random projection
    u = np.cos(y) * np.cos(x)
    v = 0* u

    px, py,_ = pressure_solve_2d(u, v, 1.0, 1.0)

    ut, vt = u-px, v-py
    d1 = divergence(ut, vt)
    np.testing.assert_almost_equal(px,np.cos(y) * np.cos(x)/2) 
    np.testing.assert_almost_equal(py,-np.sin(y) * np.sin(x)/2) 
    np.testing.assert_almost_equal(d1, 0.0)


    # random data
    u = (np.cos(y) + np.cos(29*y)) * np.cos(x) + np.cos(30*y) * np.sin(30*x)
    v = np.random.random_sample(u.shape)
    v[0] = 0
    v[-1] = 0

    px, py,_ = pressure_solve_2d(u, v, 1.0, 1.0)

    ut, vt = u-px, v-py
    d1 = divergence(ut, vt)
    np.testing.assert_almost_equal(d1, 0.0)

    # data over another interval
    scalex = 20.0
    scaley = 3
    x  = np.arange(102)/102 *2* pi * scalex
    y = np.linspace(0,pi,51) * scaley
    y, x= np.meshgrid(y, x, indexing='ij')

    u = np.cos(y/scaley) * np.cos(x/scalex)
    v = np.sin(2*y/scaley) * np.cos(x/scalex)

    px, py,_ = pressure_solve_2d(u, v, scalex=scalex, scaley=scaley)

    ut, vt = u-px, v-py
    d1 = divergence(ut, vt, scalex=scalex, scaley=scaley)
    np.testing.assert_almost_equal(d1, 0.0)




class Linear2DSolver(object):
    L = 3
    baroclinic_vars = ['u', 'v', 't']
    axes = ('y', 'x')

    params = dict(beta=1.0, f=0.0)

    def __init__(self):
        self.grid = Grid(scales=[1.0, 1.0], size=[100, 100])

    @property
    def ix(self):
        n = len(self.baroclinic_vars)
        d = {k: slice(i, n * self.L, n)
             for i, k in enumerate(self.baroclinic_vars)}

        for i, k in enumerate(self.barotropic_vars):
            d[k] = len(self.baroclinic_vars) * self.L + i

        return d

    def f(self, q):
        """X direction fluxes"""

        u = q[self.ix['u']]
        v = q[self.ix['v']]
        t = q[self.ix['t']]

        fq = np.zeros_like(q)
        fu = fq[self.ix['u']]
        fv = fq[self.ix['v']]
        ft = fq[self.ix['t']]

        # use pressure gradient for barotropic mode
        fu[0] = +t[0]

        fu[1] = -t[1]
        fu[2] = -t[2]

        ft[1] = -u[1]
        ft[2] = -u[2] / 4

    def g(self, q):
        """Y direction fluxes"""
        u = q[self.ix['u']]
        v = q[self.ix['v']]
        t = q[self.ix['t']]

        fq = np.zeros_like(q)
        fu = fq[self.ix['u']]
        fv = fq[self.ix['v']]
        ft = fq[self.ix['t']]

        fv[0] = +t[0]
        fv[1] = -t[1]
        fv[2] = -t[2]

        ft[1] = -v[1]
        ft[2] = -v[2] / 4

    def pressure_solve(self, q, dt=1.0):
        """Pressure solver

        Parameters
        ----------
        q: [neq, ny, nx]
            unghosted data
        dt: float, optional
            time step
        """
        p0 = q[self.ix['t']][0]
        u = q[self.ix['u']][0]
        v = q[self.ix['v']][0]

        px, py, p = self.grid.pressure_solve_2d(u, v)
        u[:] = u - px
        v[:] = v - py
        p0[:] = p/dt





if __name__ == '__main__':
    test_pressure_solve_2d()
