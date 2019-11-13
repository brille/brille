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
# trellis_r = 1/(np.arange(1,30002,5000)**(1/3))
max_sizes = 1/np.arange(100, 4000, 500)
# trellis_r = 1/(np.arange(1, 15002, 5000)**(1/3))
trellis_r = 2**np.arange(-1, 2.1, 0.5)

smat, tmat = np.meshgrid(max_sizes, trellis_r)

# print('pybind object creation times:')
# it = np.nditer([smat, tmat, None, None])
# for s, t, c, e in it:
#     tictoc.tic()
#     # while not (tictoc.elapsed() > 10 or tictoc.relative_uncertainty() < 0.1):
#     while not (tictoc.elapsed() > 1 or tictoc.relative_uncertainty() < 0.5):
#         brille.BZMeshQ(bz, max_size=s, lattice_ratio=t)
#         tictoc.toc()
#     c[...] = tictoc.average()
#     e[...] = tictoc.uncertainty()
#     print('size {:6.4f} ratio {:6.4f}: <time>={:5.3f}({:5.3f}) sum(time)={:5.2f} sec'.format(s,t,c,e,tictoc.elapsed()))
# creation = it.operands[2]
# creation_uncertainty = it.operands[3]
#
# pp.figure()
# pcolormesh(smat, tmat, np.log10(creation))
# pp.xlabel(r'maximum tetrahedral volume / $\AA^{-3}$')
# pp.ylabel(r'location trellis spacing / maximum tetrahedral circumsphere radius')
#
# pp.figure()
# for i, value in enumerate(trellis_r):
#     pp.errorbar(1/max_sizes, creation[i], yerr=creation_uncertainty[i], marker='o', linestyle='-', label='{:8.6f}'.format(value))
# pp.xlabel(r'inverse maximum tetrahedral volume / $\AA^3$')
# pp.ylabel(r'BZMeshQ object creation time / s')
# pp.legend(title='[trellis spacing]/[maximum circumsphere radius]')
#
# pp.figure()
# for i, value in enumerate(max_sizes):
#     pp.errorbar(1/trellis_r**3, creation[:,i], yerr=creation_uncertainty[:,i], marker='o', linestyle='-', label='{:8.6f}'.format(value))
# pp.xlabel(r'[maximum tetrahedral circumsphere radius]/[trellis spacing]')
# pp.ylabel(r'BZMeshQ object creation time / s')
# pp.legend(title=r'maximum tetrahedral volume / $\AA^{-1}$')
#
# pp.show()

# Next investigate now choice of maximum tetrahedra size and trellis size
# effect point-interpolation time
qx, qy = np.mgrid[1-0.3:1+0.3:0.005, 1-0.3:1+0.3:0.005]

# inshape = qx.shape
# px = cen2corner(qx)
# py = cen2corner(qy)

qx = qx.reshape(qx.size, 1)
qy = qy.reshape(qy.size, 1)
qxyz = np.concatenate((qx, qy, 0*qx), axis=1)

nb = load_interpolation_data('nb')

interpolation_trellis_ratios = 2**np.arange(-3,1.51,0.1)

print('Interpolation of {} points:'.format(qx.size))
it = np.nditer([interpolation_trellis_ratios, None, None])
for tr, i, ie in it:
    se = BrEu(nb, mesh=True, max_size=0.01, lattice_ratio=tr)
    tictoc.tic()
    while not (tictoc.elapsed() > 5 or tictoc.relative_uncertainty() < 0.01):
    # while not (tictoc.elapsed() > 1 or tictoc.relative_uncertainty() < 0.5):
        se(qxyz)
        tictoc.toc()
    i[...] = tictoc.average()
    ie[...] = tictoc.uncertainty()
    print('ratio {:6.4f}: <time>={:5.3f}({:5.3f}) sum(time)={:5.2f} sec'.format(tr,i,ie,tictoc.elapsed()))
interpolation_time = it.operands[1]
interpolation_uncertainty = it.operands[2]

pp.figure()
pp.errorbar(interpolation_trellis_ratios, interpolation_time, yerr=interpolation_uncertainty, marker='o', linestyle='-')
pp.xlabel(r'trellis spacing / maximum-circumsphere radius')
pp.ylabel('interpolation of {} points / s'.format(qx.size))

pp.figure()
pp.errorbar((4/3*np.pi)/interpolation_trellis_ratios**3, interpolation_time, yerr=interpolation_uncertainty, marker='o', linestyle='-')
pp.xlabel(r'$\rho_t$ / maximum-circumsphere')
pp.ylabel('interpolation of {} points / s'.format(qx.size))

pp.show()
