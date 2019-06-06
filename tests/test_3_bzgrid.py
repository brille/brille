#!/usr/bin/env python3
import os, sys, unittest
import numpy as np

# We need to tell Python where it can find the symbz module.
addpath = os.getcwd()
# It's either in the working directory where python was called or in
# a sub-directory (called Debug using Visual Studio under Windows)
if os.path.exists('Debug'):
    addpath += "\\Debug"
sys.path.append(addpath)

from importlib import util

if util.find_spec('symbz') is not None and util.find_spec('symbz._symbz') is not None:
    import symbz as s
elif util.find_spec('_symbz') is not None:
    import _symbz as s
else:
    raise Exception("symbz module not found!")

def make_drbz(a,b,c,al=np.pi/2,be=np.pi/2,ga=np.pi/2):
    d = s.Direct(a,b,c,al,be,ga)
    r = d.star()
    bz = s.BrillouinZone(r)
    return (d,r,bz)

hasmpl  = util.find_spec('matplotlib') is not None
hasmpl &= util.find_spec('mpl_toolkits') is not None
# protect against trying to load a submodule of a non-existant module
if hasmpl:
    hasmpl &= util.find_spec('mpl_toolkits.mplot3d') is not None
if hasmpl:
    import matplotlib.pyplot as pp
    from mpl_toolkits.mplot3d.axes3d import Axes3D

def plot_points(x,title=''):
    if hasmpl:
        fig = pp.figure()
        ax = Axes3D(fig)
        ax.scatter(x[:,0],x[:,1],x[:,2],s=10)
        ax.set_title(title)
        pp.show()
def plot_points_with_lines(x,y,title=''):
    if hasmpl:
        fig = pp.figure()
        ax = Axes3D(fig)
        ax.plot(y[:,0],y[:,1],y[:,2])
        ax.scatter(x[:,0],x[:,1],x[:,2],s=10)
        ax.set_title(title)
        pp.show()



class BrillouinZoneGrid (unittest.TestCase):
    def test_a_init_unit_cube(self):
        d,r,bz = make_drbz(1,1,1)
        with self.assertRaises(RuntimeError):
            s.BZGridQ(bz,4)
        with self.assertRaises(RuntimeError):
            s.BZGridQ(bz,[4,3])
        with self.assertRaises(RuntimeError):
            s.BZGridQ(bz,[[2,2,2]])
        Ntuple = (2,2,2)
        bzg = s.BZGridQ(bz, Ntuple)
        hkl = bzg.rlu
        self.assertEqual(hkl.ndim,2)
        self.assertEqual(hkl.shape[0], np.prod(2*np.array(Ntuple)+1) )
        self.assertEqual(hkl.shape[1],3)
        xyz = bzg.invA # the positions of the points in Å⁻¹
        self.assertEqual(hkl.shape,xyz.shape)
        self.assertAlmostEqual( np.abs(xyz - 2*np.pi*hkl).sum(), 0 )
    # def test_b_plot_unit_cube(self):
    #     d,r,bz = make_drbz(2*np.pi,2*np.pi,2*np.pi)
    #     step = np.array((0.2,0.2,0.2))
    #     bzg = s.BZGrid(bz, step, False)
    #     numpoints = (2*np.ceil(1/(2*step))+2).prod()
    #     plot_points( bzg.hkl )
    #     self.assertEqual( bzg.hkl.shape[0], numpoints )
    # def test_a_init_hexagonal(self):
    #     d,r,bz = make_drbz(3,3,3,np.pi/2,np.pi/2,2*np.pi/3)
    #     with self.assertRaises(RuntimeError):
    #         s.BZGrid(bz,4)
    #     with self.assertRaises(RuntimeError):
    #         s.BZGrid(bz,[4,3])
    #     with self.assertRaises(RuntimeError):
    #         s.BZGrid(bz,[[2,2,2]])
    #     Ntuple = (20,20,1)
    #     bzg = s.BZGrid(bz, Ntuple)
    def test_b_plot_hexagonal(self):
        d,r,bz = make_drbz(3,3,3,np.pi/2,np.pi/2,2*np.pi/3)
        Ntuple = (20,20,0)
        bzg = s.BZGridQ(bz, Ntuple)
        # plot_points( bzg.invA        ,'full grid')
        plot_points_with_lines( bzg.mapped_invA, bz.vertices_invA ,'mapped grid')
    # def test_c(self):
    #     d,r,bz = make_drbz(1,1,1, 2*np.pi/3, 2*np.pi/3, np.pi/3)
    #     bzg = s.BZGrid(bz, (5,5,5))
    #     plot_points_with_lines( bzg.mapped_invA, bz.vertices_invA, 'rhomb')
    def test_d_copying(self):
        d,r,bz = make_drbz(3,3,3,np.pi/2,np.pi/2,2*np.pi/3)
        dtuple = (0.2,0.2,0.3)
        bzg0 = s.BZGridQ(bz,dtuple,False)
        bzg1 = s.BZGridQ( bzg0.BrillouinZone, bzg0.halfN)
        #
        # and for 4D grids:
        fd0 = s.BZGridQE(bz, (0.,1.,10.), (2,2,2))
        fd1 = s.BZGridQE(fd0.BrillouinZone, fd0.spec, fd0.halfN)

    def test_j_data_sum(self):
        d_lat = s.Direct((3, 3, 3), np.pi/2*np.array((1, 1, 4/3)))
        r_lat = d_lat.star()
        b_zone = s.BrillouinZone(r_lat)
        bz_grid = s.BZGridQ(b_zone, halfN=(3, 3, 3))
        q_pts = bz_grid.mapped_rlu
        # fill with 1 at each mapped point
        bz_grid.fill(np.ones(q_pts.shape[0]))
        self.assertTrue(np.isclose(bz_grid.sum_data(0), q_pts.shape[0]))
        self.assertTrue(np.isclose(bz_grid.sum_data(1), np.ones(q_pts.shape[0])).all())
        # fill with |q| at each mapped point
        mod_q = np.abs(bz_grid.mapped_invA)
        bz_grid.fill(mod_q)
        self.assertTrue(np.isclose(bz_grid.sum_data(0), mod_q.sum(0)).all())
        self.assertTrue(np.isclose(bz_grid.sum_data(1), mod_q.sum(1)).all())


if __name__ == '__main__':
  unittest.main()
