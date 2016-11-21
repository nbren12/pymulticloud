"""Module for computing the mass and loading matrices for various spectral
truncations

The inner products are preformed using Gaussian quadrature

"""
from functools import reduce
from itertools import product
from numpy import sin, cos, pi, sqrt
from scipy.integrate import quadrature
from numpy.polynomial import Legendre

def fprod(*args):
    def f(x, args=args):
        return reduce(lambda a, b: a*b(x), args, 1.0)
    return f

def basis_quad_prod(*args):
    d = {}
    for ijk in product(*args):
        f = fprod(*(a[i] for a, i in zip(args, ijk)))
        res = quadrature(f, 0, pi,tol=1e-14)[0]/pi

        d[ijk] = res

    return d


def calculate_basis(L=2):
    n = L + 1
    phi = {m: lambda x, m=m: sqrt(2) * sin(x * m) for m in range(1, n)}

    phip = {m: lambda x, m=m: sqrt(2)*cos(x * m) for m in range(1, n)}
    phip[0] = lambda x: 1.0

    phipp = {m: lambda x,m=m: sqrt(2) * sin(x * m) for m in range(1, n)}

    return phi, phip, phipp


if __name__ == '__main__':
    phi, phip, phipp =calculate_basis(2)
    print(basis_quad_prod(phip,phi,phipp))
    print(basis_quad_prod(phip,phip, phip))

    
