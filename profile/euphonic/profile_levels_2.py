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

nb = load_interpolation_data('nb')

q = np.random.rand(10000,3)*10

senb = BrEu(nb, mesh=True, parallel=True, max_size=0.0001, num_levels=5)

tictoc.tic()
while not (tictoc.elapsed() > 60 or tictoc.relative_uncertainty() < 0.05):
    senb.w_q(q, interpolate=False)
    tictoc.toc()
eu_time_point = tictoc.average()/q.shape[0]*10**6
eu_delt_point = tictoc.uncertainty()/q.shape[0]*10**6

nthreads = np.arange(1,os.cpu_count()//2+1,dtype='int')
it = np.nditer([nthreads, None, None], op_dtypes=('int','double','double'))
for n, t, u in it:
    tictoc.tic()
    while not (tictoc.elapsed() > 20 or tictoc.relative_uncertainty() < 0.05):
        senb.w_q(q, threads=n)
        tictoc.toc()
    print('{} threads: {:6.3f}/{} -> {:5.3f} +/- {:5.3f}'.format(n, tictoc.elapsed(), tictoc.times_called, tictoc.average(), tictoc.uncertainty()))
    t[...]=tictoc.average()
    u[...]=tictoc.uncertainty()
time_point = it.operands[1]/q.shape[0]*10**6
uncertainty_point = it.operands[2]/q.shape[0]*10**6

pp.figure()
pp.errorbar([1],eu_time_point, yerr=eu_delt_point, marker='s', label='Euphonic')
pp.errorbar(nthreads, time_point, yerr=uncertainty_point, marker='o', color='k', linestyle='')
pp.xlabel('Number of threads')
pp.ylabel('t/μs')
pp.title(r'Time to calculate $\epsilon$(Q) and $\omega$(Q) per Q point')
pp.legend()


ticks = [1/x for x in range(1,7)]
labs = [r'$'+'{}'.format(x)+r'^{-1}$' for x in range(1,7)]
pp.figure()
pp.errorbar(1/nthreads, time_point/max(time_point), yerr=uncertainty_point/max(time_point), marker='o', color='k', linestyle='')
pp.setp(pp.gca(),'xticks',ticks,'xticklabels',labs,'yticks',ticks,'yticklabels',labs)
pp.xlabel('inverse number of threads')
pp.title('Relative interpolation time per point')

pp.show()
