#!/usr/bin/env python3
"""Run tests of the interpolation functionality."""
import os
import sys
import unittest
import numpy as np
from pathlib import Path
from importlib.util import find_spec
addpaths = [Path(), Path('..')] # we might be in .../build/wrap/
config = os.environ.get('CMAKE_CONFIG_TYPE') # set by ctest -C <cfg>
if config:
  for path in addpaths:
    if Path(path, config).exists():
      addpaths.append(Path(path,config))
sys.path[:0] = [str(path.absolute()) for path in addpaths]

if find_spec('_brille') is not None:
    import _brille as s
elif find_spec('brille') is not None and find_spec('brille._brille') is not None:
    import brille as s
else:
    abspaths = [str(path.absolute()) for path in addpaths]
    raise Exception("brille module not found in {}!".format(abspaths))

def sqwfunc_ones(Q):
    """S(Q,W) function that is all ones."""
    return 1.0+0*(Q[:, 0]+Q[:, 1]+Q[:, 2])


def sqwfunc_x(Q):
    """S(Q,W) function that returns Qx."""
    return Q[:, 0]


def sqwfunc_y(Q):
    """S(Q,W) function that returns Qy."""
    return Q[:, 1]


def sqwfunc_z(Q):
    """S(Q,W) function that returns Qz."""
    return Q[:, 2]


def sqwfunc_xy(Q):
    """S(Q,W) function that returns Qx+Qy."""
    return Q[:, 0]+Q[:, 1]


def sqwfunc_xyz(Q):
    """S(Q,W) function that returns Qx+Qy+Qz."""
    return Q[:, 0]+Q[:, 1]+Q[:, 2]


def vecfun_ident(Q):
    """S(Q,W) function that returns vector Q."""
    return Q


def vecfun_rotz(Q, θ=0):
    """S(Q,W) function that returns vector Q rotated around the z-axis."""
    x = Q[:, [0]]
    y = Q[:, [1]]
    z = Q[:, [2]]
    c = np.cos(θ)
    s = np.sin(θ)
    return np.concatenate((x*c-y*s, x*s+y*c, z), axis=1)


def vecfun_rotx(Q, θ=0):
    """S(Q,W) function that returns vector Q rotated around the x-axis."""
    x = Q[:, [0]]
    y = Q[:, [1]]
    z = Q[:, [2]]
    c = np.cos(θ)
    s = np.sin(θ)
    return np.concatenate((x, y*c-z*s, y*s+z*c), axis=1)


def matfun_ident(Q):
    """S(Q,W) function that returns a flattend matrix with Q along its diagonal."""
    sh = Q.shape
    z = np.ndarray((sh[0], sh[1]*sh[1]), dtype=Q.dtype)
    for i in range(sh[0]):
        z[i, :] = np.diag(Q[i, :]).flatten()
    return z


def modefun(Q):
    """S(Q,W) function that returns multiple modes."""
    sh = Q.shape
    modes = np.ndarray((sh[0], sh[1], sh[1], 3), dtype=Q.dtype)
    for i in range(sh[0]):
        modes[i, :, :, 0] = np.diag(Q[i, :])
        modes[i, 0, :, 1] = Q[i, :]
        modes[i, :, 0, 2] = -Q[i, :]
        # modes[i, sh[1]-1, :, 3] = 1-Q[i, :]
        # modes[i, :, sh[1]-1, 4] = 1+Q[i, :]
    return modes


def complex_scalar(Q):
    """S(Q,W) function that returns a complex scalar number."""
    return (Q[:, 0]-Q[:, 1])+(Q[:, 2]+Q[:, 1])*1j


def setup_grid(iscomplex=False, halfN=(2, 2, 2)):
    """Create a grid object for interpolating."""
    rlat = s.Reciprocal((1, 1, 1), (90, 90, 90))
    bz = s.BrillouinZone(rlat)
    max_volume = rlat.volume/(8*np.prod(halfN))
    if iscomplex:
        bzg = s.BZTrellisQcc(bz, max_volume)
    else:
        bzg = s.BZTrellisQdd(bz, max_volume)
    return bzg


def define_Q_points(rand=False, N=10):
    """Define a number of Q points for use in inerpolating."""
    if rand:
        # the tests only work if the Q points are already within the first BZ
        Q = np.random.rand(N, 3)-0.5  # (-0.5,0.5)
    else:
        Q = np.array([[0, 0, 0],
                      [0.1, 0, 0],
                      [0, 0.1, 0],
                      [0, 0, 0.1]], dtype='double')
    return Q


def fe_dispersion(Q, p=(-16, 0.01)):
    """Calculate the dispersion relationship for spinwaves in iron."""
    J = p[0]  # exchange, meV
    d = p[1]  # anisotropy
    # the dispersion relationship:
    q_dep = np.cos(np.pi*Q[:, 0])*np.cos(np.pi*Q[:, 1])*np.cos(np.pi*Q[:, 2])
    w = d + 8*J*(1-q_dep)
    return w


def fe_analytic(Q, E, p):
    """Calculate the analytic function for S(Q,ω) for spinwaves in iron."""
    g = p[2]  # mode-lifetime, meV
    T = p[3]  # Temperature, K
    A = p[4]  # scale-factor
    w = fe_dispersion(Q, p)
    S = A/np.pi*(E/1-np.exp(-11.602*E/T))*(4*g*w)/((E**2-w**2)**2+4*(g*E)**2)
    return S


class Interpolate (unittest.TestCase):
    """Class to perform unit tests related to interpolation."""

    def test_a_norm(self):
        """Test interpolation normalisation."""
        bzg = setup_grid()
        Qi = define_Q_points()
        bzg.fill(sqwfunc_ones(bzg.rlu), (1,), bzg.rlu, (0,3))
        interpolated_ones, _ = bzg.interpolate_at(Qi)
        self.assertTrue(np.allclose(interpolated_ones, 1.))

    def test_b_x(self):
        """Test with data as Qx."""
        bzg = setup_grid()
        Qi = define_Q_points()
        bzg.fill(sqwfunc_x(bzg.rlu), (1,), bzg.rlu, (0,3))
        intres, _ = bzg.interpolate_at(Qi)
        self.assertTrue(np.isclose(np.squeeze(intres), sqwfunc_x(Qi)).all())

    def test_b_y(self):
        """Test with data as Qy."""
        bzg = setup_grid()
        Qi = define_Q_points()
        bzg.fill(sqwfunc_y(bzg.rlu), (1,), bzg.rlu, (0,3))
        intres, _ = bzg.interpolate_at(Qi)
        self.assertTrue(np.isclose(np.squeeze(intres), sqwfunc_y(Qi)).all())

    def test_b_z(self):
        """Test with data as Qz."""
        bzg = setup_grid()
        Qi = define_Q_points()
        bzg.fill(sqwfunc_z(bzg.rlu), (1,), bzg.rlu, (0,3))
        intres, _ = bzg.interpolate_at(Qi)
        self.assertTrue(np.isclose(np.squeeze(intres), sqwfunc_z(Qi)).all())

    def test_c_xy(self):
        """Test with data as Qx+Qy."""
        bzg = setup_grid()
        Qi = define_Q_points()
        bzg.fill(sqwfunc_xy(bzg.rlu), (1,), bzg.rlu, (0,3))
        intres, _ = bzg.interpolate_at(Qi)
        self.assertTrue(np.isclose(np.squeeze(intres), sqwfunc_xy(Qi)).all())

    def test_c_xyz(self):
        """Test with data as Qx+Qy+Qz."""
        bzg = setup_grid()
        Qi = define_Q_points()
        bzg.fill(sqwfunc_xyz(bzg.rlu), (1,), bzg.rlu, (0,3))
        intres, _ = bzg.interpolate_at(Qi)
        self.assertTrue(np.isclose(np.squeeze(intres), sqwfunc_xyz(Qi)).all())

    def test_d_vec_ident(self):
        """Test with data as vector Q."""
        bzg = setup_grid()
        Qi = define_Q_points()
        bzg.fill(vecfun_ident(bzg.rlu), (0,3), bzg.rlu, (0,3))
        intres, _ = bzg.interpolate_at(Qi)
        self.assertTrue(np.isclose(intres, Qi).all())

    def test_d_vec_rotx(self):
        """Test with data as vector Q rotated about x."""
        bzg = setup_grid()
        Qi = define_Q_points()
        ang = np.pi/3
        bzg.fill(vecfun_rotx(bzg.rlu, ang), (0,3), bzg.rlu, (0,3))
        intres, _ = bzg.interpolate_at(Qi)
        antres = vecfun_rotx(Qi, ang)
        self.assertTrue(np.isclose(intres, antres).all())

    def test_d_vec_rotz(self):
        """Test with data as vector Q rotated about z."""
        bzg = setup_grid()
        Qi = define_Q_points()
        ang = 3*np.pi/5
        bzg.fill(vecfun_rotz(bzg.rlu, ang), (0,3), bzg.rlu, (0,3))
        intres, _ = bzg.interpolate_at(Qi)
        antres = vecfun_rotz(Qi, ang)
        self.assertTrue(np.isclose(intres, antres).all())

    def test_e_mat_ident(self):
        """Test with data as matrix Q."""
        bzg = setup_grid()
        Qi = define_Q_points()
        bzg.fill(matfun_ident(bzg.rlu), (0,0,9), bzg.rlu, (0,3))
        intres, _ = bzg.interpolate_at(Qi)
        antres = matfun_ident(Qi)
        self.assertTrue(np.isclose(intres, antres).all())

    def test_f_complex_scalar(self):
        """Test with data as complex-valued scalars."""
        bzg = setup_grid(iscomplex=True)
        Qi = define_Q_points()
        complexQ = complex_scalar(bzg.rlu)
        bzg.fill(complexQ, (1,), complexQ, (1,))
        # complex eigenvalues have fixed phase while complex eigenvectors
        # are allowed arbitrary phase.
        val, vec = bzg.interpolate_at(Qi)
        self.assertTrue(np.isclose(np.squeeze(val), complex_scalar(Qi)).all())

    # def test_g_big_grid(self):
    #     """Test with a large grid."""
    #     bzg = setup_grid(halfN=(20, 20, 20))
    #     Qi = define_Q_points(rand=True, N=100000)
    #     bzg.fill(matfun_ident(bzg.rlu))
    #     intres = bzg.interpolate_at(Qi, True, True, 10)
    #     antres = matfun_ident(Qi)
    #     self.assertTrue(np.isclose(intres, antres).all())

    def test_i_iron_self_consistency(self):
        """Test with data as iron spinwaves, but test only *at* grid points."""
        d = s.Direct((2.87, 2.87, 2.87), np.pi/2*np.array((1, 1, 1)), "Im-3m")
        r = d.star
        bz = s.BrillouinZone(r) # constructs an irreducible Bz by default
        bzg = s.BZTrellisQdc(bz, 0.125)
        Q = bzg.rlu
        bzg.fill(fe_dispersion(Q), (1,), bzg.rlu, (0,3))
        # The irreducible interpolation must be used here
        # since the Brillouin zone is an *irreducible* Brillouin zone!!
        intres, _ = bzg.ir_interpolate_at(Q, False, False)
        antres = fe_dispersion(Q)
        self.assertTrue(np.isclose(np.squeeze(intres), antres).all())

    def test_j_stored_data_references(self):
        """Ensure that the reference count for stored data is incremented"""
        from sys import getrefcount as r
        bzg = setup_grid()
        val = sqwfunc_xyz(bzg.rlu)
        vec = vecfun_rotz(bzg.rlu, np.pi/180*40)
        # Store the value and vector array reference counts before filling
        valr, vecr = r(val), r(vec)
        # fill the data into the grid -- this increases each array refcount
        bzg.fill(val, (1,), vec, (0,3))
        self.assertTrue(valr + 1 == r(val))
        self.assertTrue(vecr + 1 == r(vec))


if __name__ == '__main__':
    unittest.main()
