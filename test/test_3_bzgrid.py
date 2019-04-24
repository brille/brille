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

import symbz as s

def make_drbz(a,b,c,al=np.pi/2,be=np.pi/2,ga=np.pi/2):
    d = s.Direct(a,b,c,al,be,ga)
    r = d.star()
    bz = s.BrillouinZone(r)
    return (d,r,bz)

from importlib import util
hasmpl  = util.find_spec('matplotlib') is not None
hasmpl &= util.find_spec('mpl_toolkits') is not None
# protect against trying to load a submodule of a non-existant module
if hasmpl:
    hasmpl &= util.find_spec('mpl_toolkits.mplot3d') is not None
if hasmpl:
    import matplotlib.pyplot as pp
    from mpl_toolkits.mplot3d.axes3d import Axes3D

def plot_points(x):
    if hasmpl:
        fig = pp.figure()
        ax = Axes3D(fig)
        ax.scatter(x[:,0],x[:,1],x[:,2],s=10)
        pp.show()



class BrillouinZoneGrid (unittest.TestCase):
    def test_a_init_unit_cube(self):
        d,r,bz = make_drbz(1,1,1)
        with self.assertRaises(RuntimeError):
            s.BZGrid(bz,4)
        with self.assertRaises(RuntimeError):
            s.BZGrid(bz,[4,3])
        with self.assertRaises(RuntimeError):
            s.BZGrid(bz,[[2,2,2]])
        Ntuple = (2,2,2)
        bzg = s.BZGrid(bz, Ntuple)
        hkl = bzg.hkl
        self.assertEqual(hkl.ndim,2)
        self.assertEqual(hkl.shape[0], np.prod(2*Ntuple) )
        self.assertEqual(hkl.shape[1],3)
        xyz = bzg.xyz # the positions of the points in Å⁻¹
        self.assertEqual(hkl.shape,xyz.shape)
        self.assertAlmostEqual( np.abs(xyz - 2*np.pi*hkl).sum(), 0 )
    def test_b_plot_unit_cube(self):
        d,r,bz = make_drbz(2*np.pi,2*np.pi,2*np.pi)
        step = np.array((0.2,0.2,0.2))
        bzg = s.BZGrid(bz, step, False)
        numpoints = (2*np.ceil(1/(2*step))+2).prod()
        plot_points( bzg.hkl )
        self.assertEqual( bzg.hkl.shape[0], numpoints )


if __name__ == '__main__':
  unittest.main()
