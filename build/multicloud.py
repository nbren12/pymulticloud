"""Wrapper module for multicloud model fortran code.

This module uses cffi to wrap a fortran function with iso_c_bindings
"""
import os
from cffi import FFI

_uname_library_extensions = {'Linux': '.so', 'Darwin': '.dylib'}


def ldd_linux(libname):
    """List libraries linked in executable"""
    import re
    import sh
    ldd_out = sh.ldd(libname)

    libraries = {}
    for line in ldd_out.splitlines():
        match = re.match(r'\t(.*) => (.*) \(0x', line)
        if match:
            libraries[match.group(1)] = match.group(2)

    return libraries


def _open_library():
    """Open the multicloud library

    Loads the linked gfotran path rather than the anaconda one
    """

    library_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'libmulticloud' + _uname_library_extensions[os.uname().sysname])

    if os.uname().sysname == 'Linux':
        libs = ldd(library_path)
        ffi.dlopen(libs['libgfortran.so.3'])
    _mc_library = ffi.dlopen(library_path)
    return _mc_library


def _as_pointer(numpy_array):
    assert numpy_array.flags[
        'F_CONTIGUOUS'], "array is not contiguous in memory (Fortran order)"
    return ffi.cast("double*", numpy_array.__array_interface__['data'][0])


ffi = FFI()
ffi.cdef("""
void multicloud_wrapper(double *, double *, double *, double *,
                        double *, double *, double *, double *,
                        double *, double *,
                        int * n , double * dt, double*  dx, double*  time, double *);

void multicloud_eq(double * fceq, double * fdeq, double * fseq);
""")

_mc_library = _open_library()


def multicloud_rhs(fc, fd, fs, u1, u2, t1, t2, teb, q, hs, dt, dx, time,
                   tebst):
    """Wrapper for multicloud model right hand side"""
    n = fc.shape[0]
    _mc_library.multicloud_wrapper(
        _as_pointer(fc), _as_pointer(fd), _as_pointer(fs), _as_pointer(u1),
        _as_pointer(u2), _as_pointer(t1), _as_pointer(t2), _as_pointer(teb),
        _as_pointer(q), _as_pointer(hs), ffi.new("int *", n),
        ffi.new("double *", dt), ffi.new("double *", dx),
        ffi.new("double *", time), _as_pointer(tebst))


def equilibrium_fractions():
    """Get equilibrium cloud fractions from multicloud model
    """
    fc = ffi.new("double *", 0)
    fd = ffi.new("double *", 0)
    fs = ffi.new("double *", 0)
    _mc_library.multicloud_eq(fc, fd, fs)

    return fc[0], fd[0], fs[0]


def test_multicloud_rhs():
    import numpy as np
    n = 100

    dt = .01
    dx = .1
    time = 0

    tebst = np.zeros(n)
    fc = np.zeros(n)
    fd = np.zeros(n)
    fs = np.zeros(n)
    u1 = np.zeros(n)
    u2 = np.zeros(n)
    t1 = np.zeros(n)
    t2 = np.zeros(n)
    teb = np.zeros(n)
    q = np.zeros(n)
    hs = np.zeros(n)

    multicloud_rhs(fc, fd, fs, u1, u2, t1, t2, teb, q, hs, dt, dx, time, tebst)


if __name__ == '__main__':
    print(equilibrium_fractions())
