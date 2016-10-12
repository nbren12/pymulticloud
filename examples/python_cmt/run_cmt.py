#!/usr/bin/env python
"""Run python multicloud model

Usage:
    run_mc.py [-r <restart_file>] [--duration=<time>] [--output_interval=<interval>] [--cfl=<float>]

Options:
    -r --restart                        use restart file to initialize run
    --cfl=<float>                       cfl [default: .1]
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

        self._du = [0,0,0]


    def init_mc(self, *args, **kwargs):

        soln, dx = self._multicloud_model.init_mc(*args, extra_vars=['scmt'],
                                                  **kwargs)

        return soln, dx

    def onestep(self, soln, time, dt, *args, **kwargs):

        soln = self._multicloud_model.onestep(soln, time, dt, *args, **kwargs)
        soln = self._cmt_step(soln, time, dt)

        return soln

    def _cmt_step(self, soln, time, dt):
        """Step of cmt model"""



        u  = soln['u'] * cmt.c
        hd  = soln['hd'] * cmt.qscale
        hc  = soln['hc'] * cmt.qscale
        lmd = soln['lmd']
        lmd = 0 * lmd + .2
        scmt = soln['scmt'].astype(np.int32)

        rates, dulow, dumid = cmt.transition_rates_array(u, hd, hc, lmd)
        scmt = cmt.stochastic_integrate_array(scmt, rates, cmt.T*time, cmt.T*(time + dt))


        u = cmt.update_cmt(u, scmt, hd, dulow, dumid, dt*cmt.T)

        soln['u'] = u /cmt.c
        soln['scmt'] = scmt.astype(np.float64)

        return soln

    def __getattr__(self, name):
        return self._multicloud_model.__getattribute__(name)

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)
    import logging
    logging.basicConfig(level=logging.INFO)


    restart_file = None
    if args['--restart']:
        shutil.copy(args['<restart_file>'], 'restart.pkl')
        restart_file = args['<restart_file>']

    solver = CmtSolver()
    main(run_duration=float(args['--duration']),
         dt_out=float(args['--output_interval']),
         restart_file = restart_file, 
         cfl=float(args['--cfl']),
         solver=solver)


