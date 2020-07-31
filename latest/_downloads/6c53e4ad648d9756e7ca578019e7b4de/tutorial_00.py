import os
import numpy as np, matplotlib as mpl, matplotlib.pyplot as pp
import brille as b, brille.plotting as bp
from brilleu.brilleu import getBrillEuObj

images_dir = os.path.join(os.path.dirname(__file__), 'images')
def dirsavefig(dirname, filename, transparent=True, **kwds):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    pp.savefig(os.path.join(dirname, filename), transparent=transparent, **kwds)

def centre_to_corner(X):
    dX = np.diff(X, axis=0)/2
    X = np.concatenate((
        X[:-1]-dX,
        X[np.newaxis, -1]-dX[np.newaxis, -1],
        X[np.newaxis, -1]+dX[np.newaxis, -1]), axis=0)
    dY = np.diff(X, axis=1)/2
    X = np.concatenate((
        X[:, :-1]-dY,
        X[:, np.newaxis, -1]-dY[:, np.newaxis, -1],
        X[:, np.newaxis, -1]+dY[:, np.newaxis, -1]), axis=1)
    return X

# use_c=False ensures the Euphonic C module is *not* used.
# but *need* to use the brille C++ module, and parallel=True ensures we do so with OpenMP
scattering_lengths = {'Na': 3.63, 'Cl': 9.577} # in fm
nacl = getBrillEuObj('NaCl', max_volume=1e-5, parallel=True, sort=True, scattering_lengths=scattering_lengths)

# Define a path through reciprocal space from (000) to (123)
xi = np.linspace(0, 1, 100)
qhkl = np.vstack((xi, 2*xi, 3*xi)).T
# Interpolate the mode energies at qhkl
w_q = nacl.w_q(qhkl)

bz = nacl.grid.BrillouinZone
# plot the first Brillouin zones for (000), (111), (022), and (113) and path
fig = pp.figure(figsize=pp.figaspect(2))
ax = fig.add_axes([0,0,1,1], projection='3d')
ax.view_init(5,300)
for g in [(0,0,0), (1,1,1), (0,2,2), (1,1,3)]:
    bp.plot(bz, units='rlu', irreducible=False, origin=g, edgecolor='gray', alpha=0, axs=ax, show=mpl.is_interactive())
bp.plot(qhkl, axs=ax, show=mpl.is_interactive())
pp.setp(ax,'xlim',[-1,2],'ylim',[-1,2],'zlim',[-1,4])
pp.setp(ax,'xlabel','h','ylabel','k','zlabel','l')
bbox = fig.bbox_inches.from_bounds(0.2,0.8,4.4,5.7)
dirsavefig(images_dir, 'nacl_123_path.png', bbox_inches=bbox)

# plot the Brillouin zone, path, irreducible path, and dispersion
irqhkl, G, rot, invrot = nacl.grid.BrillouinZone.ir_moveinto(qhkl)
fig = pp.figure(figsize=pp.figaspect(0.5))
ax0 = fig.add_subplot(1, 2, 1, projection='3d')
bp.plot(bz, units='rlu', axs=ax0, Q=irqhkl, show=mpl.is_interactive())
pp.setp(ax0,'xlabel','h','ylabel','k','zlabel','l')

ax1 = fig.add_subplot(1, 2, 2)
for y in w_q.T:
	ax1.plot(xi, y)
ax1.set_xlabel(r'$\xi$ in $\mathbf{Q}=(\xi\,2\xi\,3\xi)$')
ax1.set_ylabel(r'$\omega_i(\mathbf{Q})$ / meV')

fig.subplots_adjust(left=0,right=0.99,bottom=0.1,top=1)
dirsavefig(images_dir, 'nacl_123_disp.png')

# determine and plot the inelastic neutron scattering intesity along the path
mu, nu = np.mgrid[0:1:50j, 0.5:30:60j]
ξ = mu.reshape(mu.size,1)
En = nu.reshape(nu.size,1)
Qhkl = np.concatenate((ξ, 2*ξ, 3*ξ), axis=1)

brSQE = nacl(Qhkl, En, resfun='sho', param=1., T=1.).reshape(mu.shape)
euSQE = nacl(Qhkl, En, resfun='sho', param=1., T=1., interpolate=False).reshape(mu.shape)

mu = centre_to_corner(mu)
nu = centre_to_corner(nu)
fig = pp.figure(figsize=pp.figaspect(0.5))
ax2 = fig.add_axes([0.07,0.1,0.4,0.8])
ax2.pcolormesh(mu, nu, brSQE, shading='flat')
ax2.set_xlabel(r'$\xi$ in $\mathbf{Q}=(\xi\,2\xi\,3\xi)$')
ax2.set_ylabel(r'$E$/meV')
ax2.set_title('brille')

ax3 = fig.add_axes([0.55,0.1,0.4,0.8])
ax3.pcolormesh(mu, nu, euSQE, shading='flat')
ax3.set_xlabel(r'$\xi$ in $\mathbf{Q}=(\xi\,2\xi\,3\xi)$')
ax3.set_ylabel(r'$E$/meV')
ax3.set_title('euphonic')

dirsavefig(images_dir, 'nacl_123_sqw.png')
