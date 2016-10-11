
"""Run python multicloud model

Usage:
    run_mc.py [-r <restart_file>] [--duration=<time>] [--output_interval=<interval>]

Options:
    -r --restart                        use restart file to initialize run
    -d <time> --duration=<time>         run duration [default: 100]
    -i <time> --output_interval=<time>  run duration [default: .25]
"""
import shutil
import sys

import numpy as np

sys.path.insert(0, "../../")

# this needs to be imported before python. not sure why
# import fortran.multicloud
from python.swe.multicloud import MulticloudModel, main
from python import cmt

class CmtSolver(object):
    def __init__(self):
        "docstring"
        self._multicloud_model = MulticloudModel()
        self._multicloud_model.variables.append('scmt')

        self._du = [0,0,0]

    def onestep(self, soln, time, dt, *args, **kwargs):

        soln = self._multicloud_model.onestep(soln, time, dt, *args, **kwargs)
        soln = self._cmt_step(soln, time, dt)

        return soln

    def _cmt_step(self, soln, time, dt):
        """Step of cmt model"""


        variable_idxs = self._multicloud_model.variable_idxs

        u  = soln[variable_idxs['u']] * cmt.c
        hd  = soln[variable_idxs['hd']] * cmt.qscale
        hc  = soln[variable_idxs['hc']] * cmt.qscale
        lmd = soln[variable_idxs['lmd']]
        lmd = 0 * lmd + .2
        scmt = soln[variable_idxs['scmt']].astype(np.int32)

        rates, dulow, dumid = cmt.transition_rates_array(u, hd, hc, lmd)
        scmt = cmt.stochastic_integrate_array(scmt, rates, cmt.T*time, cmt.T*(time + dt))


        # adams bashforth second order
        if np.any(u**2 > 100**2):
            i = (u**2 > 100**2).nonzero()
            raise ValueError("%i "%scmt[i])
        u += dt * cmt.T * (23/12 * self._du[0] - 4/3 * self._du[1] + 5/12*self._du[2])


        soln[variable_idxs['u']] = u /cmt.c
        soln[variable_idxs['scmt']] = scmt.astype(np.float64)

        return soln

    def __getattr__(self, name):
        return self._multicloud_model.__getattribute__(name)

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)

    restart_file = None
    if args['--restart']:
        shutil.copy(args['<restart_file>'], 'restart.pkl')
        restart_file = args['<restart_file>']

    solver = CmtSolver()
    main(run_duration=float(args['--duration']),
         dt_out=float(args['--output_interval']),
         restart_file = restart_file, 
         solver=solver)


