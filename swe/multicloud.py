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


variables = ['q' ,'teb', 'hs', 'tebst', 'fc', 'fc', 'fd', 'fs']
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


    dt = dx * .25
    return soln, dx, dt

def main():
    import logging
    logging.basicConfig(level=logging.INFO)
    soln, dx, dt = init_mc()
    t_end = 30

    for t, soln in steps(onestep, soln, dt, t_end, dx):

    np.savez("soln.npz",
             **{key: soln[variable_idxs[key]] for key in variable_idxs})

if __name__ == '__main__':
    main()
