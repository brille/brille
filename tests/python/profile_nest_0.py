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



dlat = brille.Direct((10, 10, 10.), (90.,90.,90.), 525)
bz = brille.BrillouinZone(dlat.star)

# max_sizes = 1/np.arange(100,5000,500)
max_sizes = 1/(10**np.linspace(6,12,5))**(1/3)
num_levels = np.arange(1,5).reshape(4,1)

print('pybind object creation times:')
it = np.nditer([max_sizes, num_levels, None, None, None])
for s, l, c, e, n in it:
    tictoc.tic()
    while not (tictoc.elapsed() > 30 or tictoc.relative_uncertainty() < 0.2):
        mesh = brille.BZNestQ(bz, max_volume=s, max_branchings=l)
        tictoc.toc()
    c[...] = tictoc.average()
    e[...] = tictoc.uncertainty()
    n[...] = mesh.rlu.shape[0]
    print('size {:6.4f}: <time>={:5.3f}({:5.3f}) sum(time)={:5.2f} sec'.format(s,c,e,tictoc.elapsed()))
creation = it.operands[2]
creation_uncertainty = it.operands[3]
n_vertices = it.operands[4]

pp.figure()
for l, c, e in zip(num_levels, creation, creation_uncertainty):
    pp.errorbar(max_sizes, c, yerr=e, marker='o', label=" ".join(["{:d}".format(x) for x in l]))
pp.xlabel(r'Maximum tetrahedra volume / $\AA^{-3}$')
pp.ylabel(r'BZNestQ object creation time / s')
pp.legend(title='levels')

pp.figure()
for s, c, e in zip(max_sizes, creation.T, creation_uncertainty.T):
    pp.errorbar(num_levels, c, yerr=e, marker='o', label="{:6.4f}".format(s))
pp.xlabel('Maximum tree branchings')
pp.ylabel(r'BZNestQ object creation time / s')
pp.legend(title=r'Max volume / $\AA^{-3}$')

pp.show()
