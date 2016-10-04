"""Module for timesteppers and associated utilties"""
import logging

def steps(onestep, q, dt, t_end, *args, **kwargs):
    """Iterator with fixed time step

    Parameters
    ----------
    onestep: callable(soln, t, dt, *args, **kwargs)
    """
    t = 0

    yield t, q

    while (t < t_end - 1e-10):
        dt = min(dt, t_end-t)
        q = onestep(q, t, dt, *args)
        t+=dt

        logging.info("t = {t:.2f}".format(t=t))

        yield t, q
