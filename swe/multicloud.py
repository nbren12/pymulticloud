"""Implementation of multicloud model"""
import itertools
import sys

import numpy as np

from .two_mode_swe import f as f2m
from .tadmor_1d import periodic_bc, central_scheme
from .timestepping import steps

sys.path.insert(0, "/Users/noah/workspace/multicloudstochfmk13/")
import build.multicloud as mc


L = 3


variables = ['q' ,'teb', 'hs', 'tebst', 'fc', 'fd', 'fs']
variable_idxs = dict(zip(variables, itertools.count(2*L)))
variable_idxs['u'] = slice(0,2*L, 2)
variable_idxs['t'] = slice(1,2*L, 2)

def f(q, alpha_tld=0.1,lmd_tld=0.8, q_tld=0.9):
    """Flux function for multicloud model"""

    u = q[variable_idxs['u'],...]
    T = q[variable_idxs['t'],...]
    moist = q[variable_idxs['q'],...]

    fq = np.empty_like(q)
    f2m(q, fq=fq)

    # moisture equation
    fq[variable_idxs['q'], ...] = q_tld * (u[1] + lmd_tld * u[2]) \
                                  + moist * (u[1] + alpha_tld * u[2])


    return fq

def onestep(soln, time, dt, dx):

    # hyperbolic terms
    periodic_bc(soln)
    soln[:2*L+1,...] = central_scheme(f, soln[:2*L+1,...], dx, dt)

    # multicloud model step
    mc.multicloud_rhs(soln[variable_idxs['fc']],
                      soln[variable_idxs['fd']],
                      soln[variable_idxs['fs']],
                      soln[variable_idxs['u']][1],
                      soln[variable_idxs['u']][2],
                      soln[variable_idxs['t']][1],
                      soln[variable_idxs['t']][2],
                      soln[variable_idxs['teb']],
                      soln[variable_idxs['q']],
                      soln[variable_idxs['hs']],
                      dt, dx, time,
                      soln[variable_idxs['tebst']])

    return soln


def init_mc(n=1000, dx=50/1500):
    import numpy as np
    neq = 2*L + len(variables)
    soln = np.zeros((neq, n))

    fceq, fdeq, fseq = mc.equilibrium_fractions()

    soln[variable_idxs['fc']] = fceq
    soln[variable_idxs['fd']] = fdeq
    soln[variable_idxs['fs']] = fseq

    # initialize temperature field with small random perturbation
    soln[variable_idxs['t'],...][0] = np.random.randn(n) * .01


    dt = dx * .10
    return soln, dx, dt

def unghosted(q):
    return q[:,2:-2]

def record_array_soln(soln, t):
    """Create record array for solution"""

    soln = unghosted(soln)
    n = soln.shape[-1]

    variables_dtype = [(k, np.float64, (n,)) for k in variables] \
                      + [(k, np.float64, (L, n)) for k in ['u', 't']]
    mydt = np.dtype([('time', np.float64)] + variables_dtype)
    arr = np.zeros(1, dtype=mydt)


    for k, v in variable_idxs.items():
        arr[k] = soln[v]

    arr['time'] = t

    return arr



def main():
    import logging
    logging.basicConfig(level=logging.INFO)
    soln, dx, dt = init_mc()
    t_end = 200

    dt_out = 1.0
    t_out = 0.0
    i_out = 1
    i_file = 0

    arr  = record_array_soln(soln, 0.0)

    nbuf = 10
    output = np.zeros(nbuf, dtype=arr.dtype)
    output[0] = arr

    for t, soln in steps(onestep, soln, dt, t_end, dx):

        if i_out >= nbuf:
            cur_file_name = "soln%05d.npz"%i_file
            logging.info("Saving data to file `{1}` at t={0}".format(t, cur_file_name))
            np.savez(cur_file_name, output)
            i_file += 1
            i_out = 0

        if t > t_out:
            logging.info("Storing output data at t={0}".format(t))
            output[i_out] = record_array_soln(soln, t)
            t_out += dt_out
            i_out += 1


    np.savez("soln%05d.npz"%i_file, output[:i_out])

def test_record_array_soln():

    soln, dx, dt = init_mc()
    arr = record_array_soln(soln, 0.0)

if __name__ == '__main__':
    main()
