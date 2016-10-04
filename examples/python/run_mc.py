"""Run python multicloud model

Usage:
    run_mc.py [-r <restart_file>]

Options:
    -r --restart   use restart file to initialize run
"""
import shutil
import sys

sys.path.insert(0, "../../")
from python.swe.multicloud import main

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)

    if args['--restart']:
        shutil.copy(args['<restart_file>'], 'restart.pkl')

    main()
