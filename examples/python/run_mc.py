"""Run python multicloud model

Usage:
    run_mc.py [-r <restart_file>] [--duration=<time>]

Options:
    -r --restart                 use restart file to initialize run
    -d <time> --duration=<time>  run duration [default: 100]
"""
import shutil
import sys

sys.path.insert(0, "../../")

# this needs to be imported before python. not sure why
# import fortran.multicloud
from python.swe.multicloud import main

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)

    if args['--restart']:
        shutil.copy(args['<restart_file>'], 'restart.pkl')

    main(run_duration=int(args['--duration']))
