#!/usr/bin/env python3
"""Run tests of the sort_perm functionality."""
import os
import sys
import unittest
try:
    from importlib.util import find_spec
except ImportError:
    raise Exception("This test requires importlib.util (and Python3)")
from numpy.random import permutation
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


def make_rbz(lengths, angles=np.pi/2*np.array((1, 1, 1))):
    """Make a reciprocal lattice and Brillouin zone."""
    r_lat = s.Reciprocal(lengths, angles)
    b_z = s.BrillouinZone(r_lat)
    return (r_lat, b_z)


# pylint: disable=c0103
def mode_dispersion(Q):
    """Make a crossing three-mode dispersion relationship in Q."""
    q_dep = np.cos(np.pi*Q[:, 0])*np.cos(np.pi*Q[:, 1])*np.cos(np.pi*Q[:, 2])
    accoustic = 1.0-q_dep
    optical_0 = 1.0-q_dep/2.0  # These modes cross near the zone boundary
    optical_1 = 1.7+q_dep/4.0  # and are nearly degenerate there
    modes = np.stack((accoustic, optical_0, optical_1))
    return np.transpose(modes, (1, 0))


def mode_vector(Q):
    """Make Q-varying vectors for three dispersing modes."""
    # The first mode has a vector which rotates around z, while moving along l:
    accoustic = np.array([[1, 1, 0] for z in np.pi*Q[:, 2]])/np.sqrt(2)
    # The second mode has a vector which rotates around x, while moving along k
    optical_0 = np.array([[1, -1, 0] for k in np.pi*Q[:, 1]])/np.sqrt(2)
    # The third mode has a vector which rotates around y, while moving along h
    optical_1 = np.array([[0, 0, 1] for h in np.pi*Q[:, 0]])/np.sqrt(2)
    vectors = np.stack((accoustic, optical_0, optical_1))
    return np.transpose(vectors, (1, 0, 2))


class ModeSorting(unittest.TestCase):
    """A TestCase object class to run tests of the sort_perm functionality."""

    def test_a_scalar_sorting(self):
        """Test sorting a random permutation of 0:5 for each grid point."""
        _, b_z = make_rbz((1, 1, 1))
        bz_grid = s.BZGridQ(b_z, (10, 0, 0))
        q_points = bz_grid.rlu
        n_q = q_points.shape[0]
        scalars = np.array([permutation(np.arange(5)) for _ in range(n_q)])
        # use a pre-sorted first point to make comparison easier
        scalars[0, :] = np.arange(5)
        # we need the scalars to represent different modes, so an extra
        # dimension is required
        scalars = scalars.reshape(n_q, 5, 1)
        bz_grid.fill(scalars, [1,])
        perm = bz_grid.centre_sort_perm()
        sort_res = np.array([x[y, :] for (x, y) in zip(scalars, perm)])
        all_close = np.isclose(np.diff(sort_res, axis=0), 0).all()
        self.assertTrue(all_close)

    def test_b_scalar_sorting(self):
        """Test sorting a random permutation of 1E(0:5) for each grid point."""
        _, b_z = make_rbz((1, 1, 1))
        bz_grid = s.BZGridQ(b_z, (10, 0, 0))
        q_points = bz_grid.rlu
        n_q = q_points.shape[0]
        scalars = np.array([permutation(10**np.arange(5)) for _ in range(n_q)])
        # use a pre-sorted first point to make comparison easier
        scalars[0, :] = 10**np.arange(5)
        # we need the scalars to represent different modes, so an extra
        # dimension is required
        scalars = scalars.reshape(n_q, 5, 1)
        bz_grid.fill(scalars, [1,])
        perm = bz_grid.centre_sort_perm()
        sort_res = np.array([x[y] for (x, y) in zip(scalars, perm)])
        all_close = np.isclose(np.diff(sort_res, axis=0), 0).all()
        self.assertTrue(all_close)

    def test_c_vector_sorting(self):
        """Test sorting randomly permuted vectors."""
        _, b_z = make_rbz((1, 1, 1))
        bz_grid = s.BZGridQ(b_z, (10, 0, 0))
        q_points = bz_grid.rlu
        n_q = q_points.shape[0]
        vectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1],
                            [1, 1, 0], [-1, 0, -1], [1, 1, 1]])
        rand_vecs = np.array([permutation(vectors) for _ in range(n_q)])
        # rand_vecs is (n_q, 6, 3), which is what we want
        bz_grid.fill(rand_vecs, [0,3])
        perm = bz_grid.centre_sort_perm()
        sort_vecs = np.array([x[y, :] for (x, y) in zip(rand_vecs, perm)])
        all_close = np.isclose(np.diff(sort_vecs, axis=0), 0).all()
        self.assertTrue(all_close)

    def test_d_scalar_vector(self):
        """Test sorting randomly permuted scalars and vectors."""
        _, b_z = make_rbz((1, 1, 1))
        bz_grid = s.BZGridQ(b_z, (10, 0, 0))
        q_points = bz_grid.rlu
        n_q = q_points.shape[0]
        scalars = 10**np.arange(3).reshape(3, 1)
        vectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        sca_vec = np.concatenate((scalars, vectors), axis=1)
        r_gsv = np.array([permutation(sca_vec) for _ in range(n_q)])
        bz_grid.fill(r_gsv, [1,3])
        perm = bz_grid.centre_sort_perm()
        sort_gsv = np.array([x[y, :] for x, y in zip(r_gsv, perm)])
        all_close = np.isclose(np.diff(sort_gsv, axis=0), 0).all()
        self.assertTrue(all_close)

    # # pylint: disable=r0914
    # def test_mode_crossing(self):
    #     """Test sorting randomly permuted crossing modes.
    #
    #     At the moment this test fails due to problems with sorting.
    #     """
    #     w_en = 1
    #     w_vc = 1
    #     _, b_z = make_rbz((1, 1, 1))
    #     bz_grid = s.BZGridQ(b_z, (15, 15, 15))
    #     q_points = bz_grid.rlu
    #     energies = mode_dispersion(q_points)
    #     vectors = mode_vector(q_points)
    #     n_pt = q_points.shape[0]
    #     n_br = energies.shape[1]
    #     energies = energies.reshape((n_pt, n_br, 1))
    #     en_vec = np.concatenate((energies, vectors), axis=2)
    #     # for sorting we need (n_pt, n_br, 4), which it is
    #     # make a randomly permuted version of the en_vec
    #     r_en_vec = np.array([x[permutation(range(n_br)), :] for x in en_vec])
    #     # ensure that we can compare the sorted version by fixing the first set
    #     r_en_vec[0, :] = en_vec[0, :]
    #     # put the random version into the grid
    #     # print(r_en_vec.shape)
    #     bz_grid.fill(r_en_vec, [1,3])
    #     # sort it, using that each branch has 1 energy and 3 vector elements
    #     perm = bz_grid.centre_sort_perm(scalar_cost_weight=w_en,
    #                                     eigenvector_cost_weight=w_vc)
    #     s_en_vec = np.array([x[y, :] for (x, y) in zip(r_en_vec, perm)])
    #     # before checking that the sorted result is correct, first check that
    #     # the sorting itself is stable:
    #     bz_grid.fill(s_en_vec, [1,3])
    #     new_perm = bz_grid.centre_sort_perm(scalar_cost_weight=w_en,
    #                                         eigenvector_cost_weight=w_vc)
    #     # each row of new_perm should be range(n_br) if the sorting is stable
    #     test_stability = np.array([(x == range(n_br)).all() for x in new_perm])
    #     self.assertTrue(test_stability.all())
    #     test_result = np.isclose(en_vec, s_en_vec)
    #     # if not test_result.all():
    #     #     np.set_printoptions(linewidth=2000)
    #     #     np.set_printoptions(threshold=sys.maxsize)
    #     #     for (t_r, a_v, b_v) in zip(test_result, en_vec, s_en_vec):
    #     #         if (~t_r).any():
    #     #             t_p = (~t_r).any(0)
    #     #             print(np.concatenate((a_v[:, t_p], b_v[:, t_p]), axis=1))
    #     #     print("\nThere are ",
    #     #           (~test_result).any(2).any(1).sum(),
    #     #           " grid points with swapped values\n")
    #     self.assertTrue(test_result.all())


if __name__ == '__main__':
    unittest.main()
