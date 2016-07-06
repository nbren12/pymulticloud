import numpy as np
from scipy.io import FortranFile

# datatypes for reading fortran file
rdtype = np.float64
idtype = np.int32


def read_header(name='header.bin'):
    header = FortranFile(name)

    n = header.read_record(idtype)[0]
    ntrunc = header.read_record(idtype)[0]
    neq = 2 * ntrunc + 1

    header.close()

    return {'n': n, 'ntrunc': ntrunc, 'neq': neq}


def read_int_snapshot(fr, head):
    return dict(scmt=fr.read_record(idtype)[None, ...])


def read_real_snapshot(fr, head):
    time = fr.read_record(rdtype)
    uc = fr.read_record(rdtype).reshape((head['n'], head['neq']))
    return dict(time=time,
                u=uc[None, :, :head['ntrunc']],
                th=uc[None, :, head['ntrunc']:2 * head['ntrunc']],
                q=uc[None, :, -1],
                teb=fr.read_record(rdtype)[None, ...],
                hc=fr.read_record(rdtype)[None, ...],
                hd=fr.read_record(rdtype)[None, ...],
                hs=fr.read_record(rdtype)[None, ...],
                fcls=fr.read_record(rdtype)[None, ...],
                fdls=fr.read_record(rdtype)[None, ...],
                fsls=fr.read_record(rdtype)[None, ...])


def fiter(fun):
    """Create an iterator from a snapshot reader"""

    def outfun(name, head):
        fr = FortranFile(name)

        while True:
            try:
                yield fun(fr, head)
            except TypeError:
                fr.close()
                return

    return outfun


def merge_read_iter(it):
    """merge iterator which yields dicts of numpy arrays"""

    first = next(it)

    out = {}
    for key in first:
        out[key] = []

    for item in it:
        for key in item:
            out[key].append(item[key])

    for key in out:
        out[key] = np.concatenate(out[key])

    return out


real_snapshot_iter = fiter(read_real_snapshot)
int_snapshot_iter = fiter(read_int_snapshot)

if __name__ == '__main__':
    head = read_header()
    iarr = merge_read_iter(int_snapshot_iter('int.bin', head))
    rarr = merge_read_iter(real_snapshot_iter('real.bin', head))

    from pylab import pcolormesh, show
    pcolormesh(rarr['u'][:, :, 3])
    show()
