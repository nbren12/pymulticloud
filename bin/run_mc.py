#!/usr/bin/env python
"""Run python multicloud model

Usage:
    run_mc.py [-r <restart_file>]
              [--duration=<time>]
              [--output_interval=<interval>]
              [--cfl=<float>]
              [--solver=<name>]
              [--solver-args=<args>]
              [--init-args=<args>]

Options:
    -r --restart                        use restart file to initialize run
    --cfl=<float>                       cfl [default: .1]
    -d <time> --duration=<time>         run duration [default: 100]
    -i <time> --output_interval=<time>  run duration [default: .25]
    --solver=<name>                     solver [default: dissip]
    --solver-args=<args>                solver initialization arguments
    --init-args=<args>                  initial condition arguments
"""
import sys
import os

root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, root_dir)

# this needs to be imported before python. not sure why
# import fortran.multicloud
from python.swe.multicloud import main
import python.swe.multicloud


def args_to_dict(args):
    if args is None:
        return {}
    else:
        return eval("dict("+args+")")

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)
    import logging
    logging.basicConfig(level=logging.INFO)


    restart_file = None
    if args['--restart']:
        restart_file = args['<restart_file>']


    solver_name = args['--solver']

    solver_args = args_to_dict(args['--solver-args'])

    ic_args =args_to_dict(args['--init-args'])

    if solver_name == 'dissip':
        solver = python.swe.multicloud.\
                 MulticloudModelDissipation(**solver_args)
    elif solver_name =='cmt':
        import python.cmt
        solver = python.cmt.CmtSolver()
    elif solver_name =='cmtnonlin':
        import python.cmt
        mcsolver = python.swe.multicloud.\
                   MulticloudModelNonlinear(dissipation=0.0, **solver_args)
        solver = python.cmt.CmtSolver(mcsolver=mcsolver)
    elif solver_name == 'nonlin':
        solver = python.swe.multicloud.MulticloudModelNonlinear(**solver_args)
    else:
        raise ValueError("Solver type `{}` is not available".format(solver_name))


    main(run_duration=float(args['--duration']),
         dt_out=float(args['--output_interval']),
         restart_file = restart_file,
         cfl=float(args['--cfl']),
         solver=solver,
         init_kwargs=ic_args)

