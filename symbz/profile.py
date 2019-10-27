from contextlib import contextmanager
from timeit import default_timer

@contextmanager
def timer():
    start = default_timer()
    elapsed = lambda: default_timer()-start
    yield lambda: elapsed()
    end = default_timer()
    elapsed = lambda: end-start

import numpy as np, matplotlib.pyplot as pp
import symbz
from symbz.euphonic import SymEu
from euphonic.data.interpolation import InterpolationData

dlat = symbz.Direct((5.,5.,5.), (90.,90.,90.), 525)
bz = symbz.BrillouinZone(dlat.star)

max_sizes = 10**(np.arange(-1,-3.5,-0.25))
trellis_r = np.arange(0.25, 2.01, 0.25)

smat, tmat = np.meshgrid(max_sizes, trellis_r)
creation_time = np.zeros_like(smat)
for i, (srow, trow) in enumerate(zip(smat,tmat)):
    for j, (s, t) in enumerate(zip(srow, trow)):
        n = 0;
        with timer() as tic:
            while tic() < 1.:
                symbz.BZMeshQ(bz, max_size=s, lattice_ratio=t)
                n += 1
        creation_time[i,j] = tic()/n

# pp.pcolormesh(max_sizes, trellis_r, creation_time)
pp.plot(max_sizes, creation_time[0])
