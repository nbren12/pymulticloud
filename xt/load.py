# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 17:29:39 2015

@author: noah
"""


import os
import numpy as np
from pylab import *

fields = [('u1', 'm/s', ),
          ('u2', 'm/s', ),
          ('theta1', 'kelvin', ),
          ('theta2', 'kelvin', ),
          ('theta_eb', 'kelvin', ),
          ('q', 'kelvin', ),
          ('hs', 'kelvin/day', ),
          ('hc', 'kelvin/day', ),
          ('hd', 'kelvin/day', ),
          ('fcls', '',),
          ('fdls', ''),
          ('fsls', '',)]


def read_snapshots(file):
    """Function for reading binary data from file"""
    nfields = len(fields)
    f = open(file, 'r')
    n = np.fromfile(f, np.int32, count=1)
    data = np.fromfile(f, dtype=np.float64)
    f.close()

    # Output dictionary
    d = {}

    data = data.reshape((-1, 1 + nfields * n))

    d['t'] = data[:,0]
    out = data[:,1:].reshape((-1, nfields, n))



    for i, (name, unit ) in enumerate(fields):
        d[name] = out[:,i,:]

    d['x'] = np.linspace(0, 40, n, endpoint=False)

    return d


class Data(object):
    """Object for openning output from multicloud model"""

    n = 1040
    dtype = np.float64

    def __init__(self, folder, **kw):
        self.folder = folder

        for k in kw:
            setattr(self, k, kw[k])

        self._snapshots = read_snapshots(self._getfile('snap_shots'))
        self.n = self._snapshots['x'].shape[0]
        

    def _getfile(self, name):
        return os.path.join(self.folder, name)

    def __getitem__(self, key):
        if key == 'fs':

            data = np.fromfile(self._getfile('evap'), self.dtype)
            return data.reshape((-1, self.n))

        else:
            return self._snapshots[key]

    def xtargs(self, key):
        return self['x'], self['t'], self[key]

    def to_xray(self):
        import xray

        das = {}
        for varname, unit in fields:
            x, t, val = self.xtargs(varname)
            das[varname] = xray.DataArray(val,
                                          coords=(t, x),
                                          dims=('time', 'x'))
            
        return xray.Dataset(das)
        
    def __setitem__(self, key, val):
        self._snapshots[key] = val

    def calculate_diagnostics(self):
        from gnl.util import fftdiff
        self['P']  = self['hd'] + .4 * self['hs']

        self['w1']  = -fftdiff(self['u1'])
        self['w2']  = -fftdiff(self['u2'])

        moisture_flux = (self['u1'] + 0.1*self['u2']) * self['q'] + \
                        (self['u1'] + 0.8*self['u2']) * 0.9

        self['dwq']  = -fftdiff(moisture_flux)
        self['GMS']  = self['dwq'] /self['P']

    @classmethod
    def xray(cls, *args, **kwargs):
        return cls(*args, **kwargs).to_xray()


def load_ident(ident, data='/scratch/noah/job.output'):
    """Get path to data from unique ident"""
    from glob import glob
    from os.path import join

    return glob(join(data, '*' + ident + '*'))[0]

def test_data():
    """run from within a data directory"""

    d = Data('/scratch/noah/job.output/tta/')

    figure()
    pcolormesh(*d.xtargs('fs'), cmap='bwr')
    colorbar()
    axis('tight')
    ylim([200,300])

    show()


def main():
    print("Saving to netcdf file")
    Data("./").to_xray().to_netcdf("data.nc")

if __name__ == '__main__':
    main()
