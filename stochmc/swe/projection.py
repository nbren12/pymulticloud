"""Module for transforming between Eulerian and spectral in the vertical

"""
import numpy as np
from numpy import pi


class VerticalGrid(object):
    """ Class for converting from spatial to spectral coordinates

    Parameters
    ----------
    num_modes : int
        The desired number of vertical modes including the barotropic mode

    Attributes
    ----------
    c2z : (n,n)
        transition matrix from spectral to spatial
    nz: int
        number of vertical levels dealiased
    vertical_grid: (nz,)
        array of vertical grid locations
    """

    def __init__(self, num_modes):
        self.num_modes = num_modes
        self._ht = pi

    @property
    def nz(self):
        """Dealiased vertical grid size"""

        return (self.num_modes + 1) * 3 // 2

    @property
    def vertical_grid(self):
        nz = self.nz
        return np.arange(nz) * pi / nz

    @property
    def c2z(self):
        z = self.vertical_grid

        return np.vstack(np.cos(m * z) for m in range(z.shape[0])).T

    @property
    def z2c(self):
        T = self.c2z
        return np.linalg.inv(T)

    @property
    def s2z(self):
        z = self.vertical_grid

        return np.vstack(np.sin(m * z) for m in range(z.shape[0])).T

    @property
    def z2s(self):
        T = self.s2z
        iT = np.linalg.inv(T[1:, 1:])

        return np.pad(iT, ((1, 0), (1, 0)), 'constant')

    def apply_vertical_transform(self, T, u, axis):
        return np.rollaxis(np.tensordot(T, u, (1, axis)), 0, start=axis+1)


    def spatial2spectral(self, u, kind, axis=0):
        """ Transform spatial field to spectral coordinates

        Parameters
        ----------
        u :
            spatial field
        kind: str
            'c' or 's' for cosine or sine variable
        """

        if kind =='c':
            T = self.z2c
        elif kind =='s':
            T = self.z2s
        else:
            raise ValueError('kind "' + kind + '" is not valid.')

        return self.apply_vertical_transform(T, u, axis)


    def spectral2spatial(self, u, kind, axis=0):
        """ Transform spatial field to spectral coordinates

        Parameters
        ----------
        u :
            spectral field
        kind: str
            'c' or 's' for cosine or sine variable
        """

        if kind =='c':
            T = self.c2z
        elif kind =='s':
            T = self.c2z
        else:
            raise ValueError('kind "' + kind + '" is not valid.')

        return self.apply_vertical_transform(T, u, axis)




    def dealias_spectral(self, u, axis=0):
        """Dealias in spectral space"""
        fu = np.rollaxis(u, axis)
        assert fu.shape[0] == self.nz
        fu[self.num_modes:, ...] = 0.0
        return np.rollaxis(fu, 0, start=axis + 1)

    def dealias_spatial(self, u, kind, axis=0):
        """ Dealias in spatial space

        Parameters
        ----------
        u :
            spatial field
        kind: str
            'c' or 's' for cosine or sine variable

        """

        fu = np.rollaxis(u, axis)
        assert fu.shape[0] == self.nz

        fu = self.spatial2spectral(fu, kind, axis=0)
        fu[self.num_modes:, ...] = 0.0
        fu = self.spectral2spatial(fu, kind, axis=0)

        return np.rollaxis(fu, 0, start=axis + 1)



def test_VerticalGrid():
    grid = VerticalGrid(2)

    u = np.random.rand(100, grid.nz, 200)

    fu = grid.spatial2spectral(u, 'c', axis=1)
    assert fu.shape == u.shape

    fu = grid.spectral2spatial(fu, 'c', axis=1)
    np.testing.assert_allclose(u, fu)

    # different shape
    u = np.random.rand(100, grid.nz)

    fu = grid.spatial2spectral(u, 'c', axis=1)
    assert fu.shape == u.shape

    fu = grid.spectral2spatial(fu, 'c', axis=1)
    np.testing.assert_allclose(u, fu)

    return 1


if __name__ == '__main__':
    test_VerticalGrid()
