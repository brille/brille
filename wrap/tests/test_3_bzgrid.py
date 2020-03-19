#!/usr/bin/env python3
"""Run tests of the bzgrid functionality."""
import os
import sys
import unittest
try:
    from importlib.util import find_spec
except ImportError:
    raise Exception("This test requires importlib.util (and Python3)")
import numpy as np

# Try to find the brille module:
# It might be in the current directory, or a sub-directory called, e.g., Debug
ADDPATH = os.getcwd()
if os.path.exists('Debug'):
    ADDPATH += "\\Debug"
sys.path.append(ADDPATH)
# Now the actual search for the module
if find_spec('brille') is not None and find_spec('brille._brille') is not None:
    import brille as s
elif find_spec('_brille') is not None:
    # pylint: disable=e0401
    import _brille as s
else:
    raise Exception("Required brille module not found!")


# pylint: disable=c0103, r0913
def make_drbz(a, b, c, al=np.pi/2, be=np.pi/2, ga=np.pi/2):
    """Make a Direct Reicprocal and BrillouinZone object."""
    d = s.Direct(a, b, c, al, be, ga)
    r = d.star
    bz = s.BrillouinZone(r)
    return (d, r, bz)


hasmpl = find_spec('matplotlib') is not None
hasmpl &= find_spec('mpl_toolkits') is not None
# protect against trying to load a submodule of a non-existant module
if hasmpl:
    hasmpl &= find_spec('mpl_toolkits.mplot3d') is not None
if hasmpl:
    import matplotlib.pyplot as pp
    from mpl_toolkits.mplot3d.axes3d import Axes3D


def plot_points(x, title=''):
    """Plot points with a title."""
    if hasmpl:
        fig = pp.figure()
        ax = Axes3D(fig)
        ax.scatter(x[:, 0], x[:, 1], x[:, 2], s=10)
        ax.set_title(title)
        pp.show()


def plot_points_with_lines(x, y, title=''):
    """Plot points with connecting lines and a title."""
    if hasmpl:
        fig = pp.figure()
        ax = Axes3D(fig)
        ax.plot(y[:, 0], y[:, 1], y[:, 2])
        ax.scatter(x[:, 0], x[:, 1], x[:, 2], s=10)
        ax.set_title(title)
        pp.show()


class BrillouinZoneGrid(unittest.TestCase):
    """The tests for the BrilloiuinZoneGrid objects."""

    def test_a_init_unit_cube(self):
        """Test a: make a unit cube."""
        _, _, bz = make_drbz(1, 1, 1)
        with self.assertRaises(RuntimeError):
            s.BZGridQ(bz, 4)
        with self.assertRaises(RuntimeError):
            s.BZGridQ(bz, [4, 3])
        with self.assertRaises(RuntimeError):
            s.BZGridQ(bz, [[2, 2, 2]])
        Ntuple = (2, 2, 2)
        bzg = s.BZGridQ(bz, Ntuple)
        hkl = bzg.grid_rlu
        self.assertEqual(hkl.ndim, 2)
        self.assertEqual(hkl.shape[0], np.prod(2*np.array(Ntuple)+1))
        self.assertEqual(hkl.shape[1], 3)
        xyz = bzg.grid_invA  # the positions of the points in Å⁻¹
        self.assertEqual(hkl.shape, xyz.shape)
        self.assertAlmostEqual(np.abs(xyz - 2*np.pi*hkl).sum(), 0)
    # def test_b_plot_unit_cube(self):
    #     d, r, bz = make_drbz(2*np.pi, 2*np.pi, 2*np.pi)
    #     step = np.array((0.2, 0.2, 0.2))
    #     bzg = s.BZGrid(bz, step, False)
    #     numpoints = (2*np.ceil(1/(2*step))+2).prod()
    #     plot_points(bzg.hkl)
    #     self.assertEqual(bzg.hkl.shape[0], numpoints)
    # def test_a_init_hexagonal(self):
    #     d, r, bz = make_drbz(3, 3, 3, np.pi/2, np.pi/2, 2*np.pi/3)
    #     with self.assertRaises(RuntimeError):
    #         s.BZGrid(bz, 4)
    #     with self.assertRaises(RuntimeError):
    #         s.BZGrid(bz,[4, 3])
    #     with self.assertRaises(RuntimeError):
    #         s.BZGrid(bz,[[2, 2, 2]])
    #     Ntuple = (20, 20, 1)
    #     bzg = s.BZGrid(bz, Ntuple)
    #
    # def test_b_plot_hexagonal(self):
    #     """Test b: plotting a hexagon."""
    #     _, _, bz = make_drbz(3, 3, 3, np.pi/2, np.pi/2, 2*np.pi/3)
    #     Ntuple = (20, 20, 0)
    #     bzg = s.BZGridQ(bz, Ntuple)
    #     plot_points_with_lines(bzg.invA, bz.vertices_invA, 'mapped')

    # def test_d_copying(self):
    #     """Test d: copying BZGrid objects."""
    #     _, _, bz = make_drbz(3, 3, 3, np.pi/2, np.pi/2, 2*np.pi/3)
    #     dtuple = (0.2, 0.2, 0.3)
    #     bzg0 = s.BZGridQ(bz, dtuple, False)
    #     bzg1 = s.BZGridQ(bzg0.BrillouinZone, bzg0.halfN)
    #     #
    #     # and for 4D grids:
    #     fd0 = s.BZGridQE(bz, (0., 1., 10.), (2, 2, 2))
    #     fd1 = s.BZGridQE(fd0.BrillouinZone, fd0.spec, fd0.halfN)

    # def test_j_data_sum(self):
    #     """Test j: sum over the contained data."""
    #     d_lat = s.Direct((3, 3, 3), np.pi/2*np.array((1, 1, 4/3)))
    #     r_lat = d_lat.star
    #     b_zone = s.BrillouinZone(r_lat)
    #     bz_grid = s.BZGridQ(b_zone, halfN=(3, 3, 3))
    #     q_pts = bz_grid.rlu
    #     # fill with 1 at each mapped point
    #     bz_grid.fill(np.ones(q_pts.shape[0]),[1,])
    #     self.assertTrue(np.isclose(bz_grid.sum_data(0), q_pts.shape[0]))
    #     qones = np.ones(q_pts.shape[0])
    #     self.assertTrue(np.isclose(bz_grid.sum_data(1), qones).all())
    #     # fill with |q| at each mapped point
    #     mod_q = np.abs(bz_grid.invA)
    #     bz_grid.fill(mod_q, [1,])
    #     self.assertTrue(np.isclose(bz_grid.sum_data(0), mod_q.sum(0)).all())
    #     self.assertTrue(np.isclose(bz_grid.sum_data(1), mod_q.sum(1)).all())


if __name__ == '__main__':
    unittest.main()
