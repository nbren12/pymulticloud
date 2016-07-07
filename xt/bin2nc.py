#!/usr/bin/env python
"""
Usage:
  bin2nc.py plot FILE FIELD
  bin2nc.py report FILE

Arguments:
  FILE                output file
  FIELD               field to plot

"""
import struct
from docopt import docopt
import numpy as np


def read_header(name='real.bin'):

    fmt = 'iidddd'
    with open(name, "rb") as f:
        s = struct.unpack_from(fmt, f.read(struct.calcsize(fmt)))

    out = dict(zip(('n', 'ntrunc', 'dx', 'T', 'L', 'abar'), s))
    out['offset'] = struct.calcsize(fmt)

    return out


def my_dtype(head):

    n = head['n']
    ntrunc = head['ntrunc']

    mydt = np.dtype([('time', np.float64),
                     ('u', np.float64, (n, ntrunc)),
                     ('th', np.float64, (n, ntrunc)),
                     ('q', np.float64, (n, )),
                     ('teb', np.float64, (n, )),
                     ('hc', np.float64, (n, )),
                     ('hd', np.float64, (n, )),
                     ('hs', np.float64, (n, )),
                     ('fcls', np.float64, (n, )),
                     ('fdls', np.float64, (n, )),
                     ('fsls', np.float64, (n, )),
                     ('scmt', np.int32, (n, ))])

    return mydt


def read_output(name):
    head = read_header(name)
    data = np.memmap(name, dtype=my_dtype(head), offset=head['offset'])

    return head, data


def report(name):
    hr = 3600.0
    day = hr * 24
    km = 1000

    head, data = read_output(name)

    vars = [field for field in data.dtype.fields if data[field].ndim > 1]

    x = np.arange(head['n']) * head['dx']/ 1000
    print(x)


if __name__ == '__main__':

    args = docopt(__doc__)

    if args['plot']:
        head, data = read_output(args['FILE'])
        from pylab import pcolormesh, show, colorbar, axis
        pcolormesh(data[args['FIELD']])
        axis('tight')
        colorbar()
        show()

    if args['report']:
        report(args['FILE'])

