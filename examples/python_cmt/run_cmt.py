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
from python.swe.multicloud import  main
from python.cmt import CmtSolver



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


