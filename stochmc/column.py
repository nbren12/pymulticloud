#!/usr/bin/env python
# coding: utf-8
"""
Run the stochastic multicloud model via the python wrapped multicloud_mod.f90

Parameters of the mc model can be changed via input.nml
"""
import matplotlib.pyplot as plt
import numpy as np
from .wrapper import multicloud
from collections import defaultdict


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

    def onestep(self,
                d,
                time=0.0,
                dt=60 / (8.3 * 3600),
                output_diags=True,
                dx=0):

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

    def run(self, d=None, nstep=100000, ioskip=100, dt=60 / 8 / 3600):
        """Run the stochastic multicloud model as a column model

        Parameters
        ----------
        d: dict, optional
            initial condition. 0 is default.
        nstep: int, optional
            number of time steps
        ioskip: int, optional
            output storage interval. Default 100.
        dt: float, optional
            time step size in nondimensional units. Default is 60 seconds.

        Returns
        -------
        dict
            Dictionary of stored prognostic and diagnostic variables

        """
        if d is None:
            d = dict(zip(self.prog_vars, self.prog_vals))

        d_output = defaultdict(list)

        for i in range(nstep):
            d, diags = self.onestep(d, dt=dt)
            if i % ioskip == 0:

                out = {}
                out.update(d)
                out.update(diags)

                for key in out:
                    d_output[key].append(out[key])

                d_output['time'].append(i * dt)

        return d_output

    def run_dataframe(self, *args, **kwargs):
        """Run multicloud model and return pandas dataframe

        See Also
        --------
        run
        """
        import pandas as pd
        d_output = self.run(*args, **kwargs)

        return pd.DataFrame.from_dict(d_output).set_index('time')


def summary_plot(df):
    """Plot summary of column model output"""
    from gnl.plots import plotiter
    pt = plotiter(
        [('fracs', ['fc', 'fd', 'fs']), ('heat', ['hc', 'hd', 'hs']),
         ('temp', ['t1', 't2', 'teb', 'q']), ('hd and q', ['hd', 'hs', 'q'])],
        yield_axis=True,
        ncol=1,
        aspect=.2,
        w=10)

    for (plotname, fields), ax in pt:
        df[fields].plot(ax=ax)


def main():
    import pandas as pd
    from gnl.plots import plotiter

    col = ColumnMulticloudModel()
    df = col.run_dataframe(nstep=100000)

    summary_plot(df)
    plt.savefig("out.pdf")


if __name__ == '__main__':
    main()
