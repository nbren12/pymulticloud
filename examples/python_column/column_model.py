
# coding: utf-8

# In[43]:

import numpy as np
from fortran import multicloud
from collections import defaultdict


# need two spaces and the end of this or it will fail

# In[44]:


class ColumnMulticloudModel(object):
    """Object for running column model simulations of the multicloud model"""
    prog_vars = ['fc', 'fd', 'fs', 'u1', 'u2', 't1', 't2', 'teb', 'q', 'hs']
    diag_vars = ['tebst', 'hc', 'hd', 'moiststab']
    
    @property
    def prog_vals(self):
        return [0.0 for v in self.prog_vars]

    @property
    def feq(self):
        return multicloud.equilibrium_fractions()
    
    def onestep(self, d, time=0.0, dt=60/(8.3*3600), output_diags=True, dx=0):

        prog_vals = [np.array([float(d[v])]) for v in self.prog_vars]
        diag_vals = [np.array([0.0]) for v in self.diag_vars]

        args = prog_vals + [dt, dx, time] + diag_vals
        multicloud.multicloud_rhs(*args)

        prog_vals = [v[0] for v in prog_vals]
        diag_vals = [v[0] for v in diag_vals]

        out = dict(zip(self.prog_vars, prog_vals))



        if output_diags:
            diags_out = dict(zip(self.diag_vars, diag_vals))
            return out, diags_out
        else:
            return out

    def run(self, d=None, nstep=100000, ioskip=100, dt=60/8/3600):
        if d is None:
            d = dict(zip(self.prog_vars, self.prog_vals))
            
        d_output = defaultdict(list)

        for i in range(nstep):
            d, diags = self.onestep(d, dt=dt)
            if i%ioskip == 0:

                out = {}
                out.update(d)
                out.update(diags)

                for key in out:
                    d_output[key].append(out[key])

                d_output['time'].append(i*dt)
                
        return d_output
                
    def run_dataframe(self, *args, **kwargs):
        d_output = self.run(*args, **kwargs)
        
        return pd.DataFrame.from_dict(d_output).set_index('time')


# In[45]:

get_ipython().magic('matplotlib')
import pandas as pd
from gnl.plots import plotiter


# In[46]:

col  = ColumnMulticloudModel()
df = col.run_dataframe(nstep=100000)


# In[47]:

pt = plotiter([
        ('fracs', ['fc', 'fd', 'fs']),
        ('heat', ['hc', 'hd', 'hs']),
        ('temp', ['t1', 't2', 'teb', 'q']),
        ('hd and q', ['hd', 'hs', 'q'])
    ], yield_axis=True, ncol=1, aspect=.2, w=10)

for (plotname, fields), ax in pt:
    df[fields].plot(ax=ax)
import matplotlib.pyplot as plt
plt.savefig("out.pdf")
