from numba import jit


def periodic_bc(u, g=2, axes=(1, )):
    """periodic bc in arbitrary dimensions
    """

    for i in axes:
        idx_in = [slice(None)] * u.ndim
        idx_out = [slice(None)] * u.ndim

        idx_in[i] = slice(g, 2 * g)
        idx_out[i] = slice(-g, None)

        u[idx_out] = u[idx_in]

        idx_in[i] = slice(-2 * g, -g)
        idx_out[i] = slice(0, g)
        u[idx_out] = u[idx_in]

@jit(nopython=True)
def fmax(a, b):
    if a > b : return a
    else: return b


@jit(nopython=True)
def fmin(a, b):
    if a < b : return a
    else: return b

@jit(nopython=True)
def minmod(a, b, c):
    mmax = fmax(fmax(a,b), c)
    mmin = fmin(fmin(a,b), c)

    if mmin > 0.0:
        return mmin
    elif mmax < 0.0:
        return mmax
    else:
        return 0.0
