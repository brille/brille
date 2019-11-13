import os
import sys
try:
    from importlib.util import find_spec
except ImportError:
    raise Exception("This test requires importlib.util (and Python3)")
from timeit import default_timer
import numpy as np, matplotlib.pyplot as pp

for module in [('euphonic', 'data', 'interpolation'), ('brille', 'euphonic')]:
    for i in range(len(module)):
        check = '.'.join(np.array(module)[0:i+1])
        if find_spec(check) is None:
            raise Exception('Required module {} not found'.format(check))
from euphonic.data.interpolation import InterpolationData
from brille.euphonic import BrEu
import brille

def load_interpolation_data(named):
    """Load a data file from the repository tests folder

    Note that this only works if the package is installed via
        python setup.py develop [--user]
    which builds the C++ library and installs it inside of the repository
    structure with simlinks in the usual package location(s).
    """
    test_spec = find_spec('brille')
    brilleroot = test_spec.submodule_search_locations[0]
    seed = os.path.join(brilleroot,'..','tests',named)
    return InterpolationData(seed)

class timer:
    def __init__(self):
        self.tic()
    def tic(self):
        self.times_called = 0
        self.start_time = default_timer()
        self.elapsed_time = 0
    def toc(self):
        self.times_called += 1
        self.elapsed_time =  default_timer()-self.start_time
        return self.elapsed_time
    def elapsed(self):
        return self.elapsed_time
    def average(self):
        return self.elapsed_time/self.times_called
    def uncertainty(self):
        return np.sqrt(self.average())/self.times_called
    def relative_uncertainty(self):
        if self.times_called < 1:
            return np.inf
        return 1/np.sqrt(self.times_called*self.elapsed_time)
tictoc = timer()

def cen2corner(Z):
    dX = np.diff(Z, axis=0)/2
    Z = np.concatenate((
        Z[:-1]-dX,
        Z[np.newaxis, -1]-dX[np.newaxis, -1],
        Z[np.newaxis, -1]+dX[np.newaxis, -1]), axis=0)
    dY = np.diff(Z, axis=1)/2
    Z = np.concatenate((
        Z[:, :-1]-dY,
        Z[:, np.newaxis, -1]-dY[:, np.newaxis, -1],
        Z[:, np.newaxis, -1]+dY[:, np.newaxis, -1]), axis=1)
    return Z

def pcolormesh(X,Y,Z, *args,**kwargs):
    return pp.pcolormesh(cen2corner(X), cen2corner(Y), Z, *args, **kwargs)



dlat = brille.Direct((5.,5.,5.), (90.,90.,90.), 525)
bz = brille.BrillouinZone(dlat.star)

# max_sizes = 1/np.arange(100,5000,500)
max_sizes = 2**np.arange(-12,-2.01,0.5)
max_sizes = 1/(10**np.linspace(6,12,5))**(1/3)

print('pybind object creation times:')
it = np.nditer([max_sizes, None, None, None])
for s, c, e, n in it:
    tictoc.tic()
    while not (tictoc.elapsed() > 10 or tictoc.relative_uncertainty() < 0.01):
        mesh = brille.BZMeshQ(bz, max_size=s, lattice_ratio=2.)
        tictoc.toc()
    c[...] = tictoc.average()
    e[...] = tictoc.uncertainty()
    n[...] = mesh.rlu.shape[0]
    print('size {:6.4f}: <time>={:5.3f}({:5.3f}) sum(time)={:5.2f} sec'.format(s,c,e,tictoc.elapsed()))
creation = it.operands[1]
creation_uncertainty = it.operands[2]
n_vertices = it.operands[3]

pp.figure()
pp.errorbar(n_vertices**2, 1000*creation, yerr=1000*creation_uncertainty, marker='o')
pp.xlabel(r'$N_v^2$')
pp.ylabel(r'BZMeshQ object creation time / ms')

pp.figure()
pp.errorbar(n_vertices, 1000*creation, yerr=1000*creation_uncertainty, marker='o')
pp.xlabel(r'$N_v$')
pp.ylabel(r'BZMeshQ object creation time / ms')

pp.show()
