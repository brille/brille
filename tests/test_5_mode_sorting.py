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

# Try to find the symbz module:
# It might be in the current directory, or a sub-directory called, e.g., Debug
ADDPATH = os.getcwd()
if os.path.exists('Debug'):
    ADDPATH += "\\Debug"
sys.path.append(ADDPATH)
# Now the actual search for the module
if find_spec('symbz') is not None and find_spec('symbz._symbz') is not None:
    import symbz as s
elif find_spec('_symbz') is not None:
    import _symbz as s
else:
    raise Exception("Required symbz module not found!")


def make_rbz(lengths, angles=np.pi/2*np.array((1, 1, 1))):
    """Make a reciprocal lattice and Brillouin zone."""
    r_lat = s.Reciprocal(lengths, angles)
    b_z = s.BrillouinZone(r_lat)
    return (r_lat, b_z)


def mode_dispersion(Q):
    """Make a crossing three-mode dispersion relationship in Q."""
    q_dep = np.cos(np.pi*Q[:, 0])*np.cos(np.pi*Q[:, 1])*np.cos(np.pi*Q[:, 2])
    accoustic = 1.0-q_dep
    optical_0 = 1.0-q_dep/2.0  # These modes cross near the zone boundary
    optical_1 = 1.7+q_dep/4.0  # and are nearly degenerate there
    optical_2 = 2.0-q_dep/8.0
    modes = np.stack((accoustic, optical_0, optical_1, optical_2))
    return np.transpose(modes, (1, 0))


def mode_vector(Q):
    """Make Q-varying vectors for three dispersing modes."""
    # The first mode has a vector which rotates around z, while moving along l:
    accoustic = np.array([[np.cos(z)+np.sin(z),
                           np.cos(z)-np.sin(z),
                           0] for z in np.pi*Q[:, 2]])
    # The second mode has a vector which rotates around x, while moving along k
    optical_0 = np.array([[0,
                           np.cos(k)+np.sin(k),
                           np.cos(k)-np.sin(k)] for k in np.pi*Q[:, 1]])
    # The third mode has a vector which rotates around y, while moving along h
    optical_1 = np.array([[np.cos(h)-np.sin(h),
                           0,
                           np.cos(h)+np.sin(h)] for h in np.pi*Q[:, 0]])
    # The fourth mode rotates around z while moving in (110)
    optical_2 = np.array([[np.cos(hk)-np.sin(hk),
                           np.sin(hk)+np.cos(hk),
                           1] for hk in np.pi*(Q[:, 0]+Q[:, 1])])
    vectors = np.stack((accoustic/np.sqrt(2),
                        optical_0/np.sqrt(2),
                        optical_1/np.sqrt(2),
                        optical_2/np.sqrt(3)))
    return np.transpose(vectors, (1, 0, 2))


class ModeSorting(unittest.TestCase):
    """A TestCase object class to run tests of the sort_perm functionality."""

    def test_a_scalar_sorting(self):
        """Test sorting a random permutation of 0:5 for each grid point."""
        _, b_z = make_rbz((1, 1, 1))
        bz_grid = s.BZGridQ(b_z, (10, 0, 0))
        q_points = bz_grid.mapped_rlu
        n_q = q_points.shape[0]
        scalars = np.array([permutation(np.arange(5)) for _ in range(n_q)])
        # use a pre-sorted first point to make comparison easier
        scalars[0, :] = np.arange(5)
        bz_grid.fill(scalars)
        perm = bz_grid.sort_perm(5, 0, 0)
        sort_res = np.array([x[y] for (x, y) in zip(scalars, perm)])
        self.assertTrue(np.isclose(sort_res, np.sort(scalars, axis=1)).all())

    def test_b_scalar_sorting(self):
        """Test sorting a random permutation of 1E(0:5) for each grid point."""
        _, b_z = make_rbz((1, 1, 1))
        bz_grid = s.BZGridQ(b_z, (10, 0, 0))
        q_points = bz_grid.mapped_rlu
        n_q = q_points.shape[0]
        scalars = np.array([permutation(10**np.arange(5)) for _ in range(n_q)])
        # use a pre-sorted first point to make comparison easier
        scalars[0, :] = 10**np.arange(5)
        bz_grid.fill(scalars)
        perm = bz_grid.sort_perm(5, 0, 0)
        sort_res = np.array([x[y] for (x, y) in zip(scalars, perm)])
        self.assertTrue(np.isclose(sort_res, np.sort(scalars, axis=1)).all())

    def test_c_vector_sorting(self):
        """Test sorting randomly permuted vectors."""
        _, b_z = make_rbz((1, 1, 1))
        bz_grid = s.BZGridQ(b_z, (10, 0, 0))
        q_points = bz_grid.mapped_rlu
        n_q = q_points.shape[0]
        vectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1],
                            [1, 1, 0], [-1, 0, -1], [1, 1, 1]])
        rand_vecs = np.array([permutation(vectors) for _ in range(n_q)])
        rand_vecs[0, :] = vectors
        # rand_vecs is (n_q, 6, 3) but we need it to be (n_q, 3, 6)
        rand_vecs = np.transpose(rand_vecs, (0, 2, 1))
        bz_grid.fill(rand_vecs)
        perm = bz_grid.sort_perm(0, 3, 0)
        sort_vecs = np.array([x[:, y] for (x, y) in zip(rand_vecs, perm)])
        # permute the sorted vectors back for comparison with vectors
        sort_vecs = np.transpose(sort_vecs, (0, 2, 1))
        all_close = np.array([np.isclose(x, vectors) for x in sort_vecs])
        self.assertTrue(all_close.all())

    def test_mode_crossing(self):
        """Test sorting randomly permuted crossing modes."""
        w_en = 1
        w_vc = 1
        _, b_z = make_rbz((1, 1, 1))
        bz_grid = s.BZGridQ(b_z, (15, 15, 15))
        q_points = bz_grid.mapped_rlu
        energies = mode_dispersion(q_points)
        vectors = mode_vector(q_points)
        n_pt = q_points.shape[0]
        n_br = energies.shape[1]
        energies = energies.reshape((n_pt, n_br, 1))
        en_vec = np.concatenate((energies, vectors), axis=2)
        # for sorting we need (n_pt, 4, n_br) instead of (n_pt, n_br, 4)
        en_vec = np.transpose(en_vec, (0, 2, 1))
        # make a randomly permuted version of the en_vec
        r_en_vec = np.array([x[:, permutation(range(n_br))] for x in en_vec])
        # ensure that we can compare the sorted version by fixing the first set
        r_en_vec[0, :] = en_vec[0, :]
        # put the random version into the grid
        bz_grid.fill(r_en_vec)
        # sort it, using that each branch has 1 energy and 3 vector elements
        perm = bz_grid.new_sort_perm(1, 3, 0, w_en, w_vc)
        s_en_vec = np.array([x[:, y] for (x, y) in zip(r_en_vec, perm)])
        # before checking that the sorted result is correct, first check that
        # the sorting itself is stable:
        bz_grid.fill(s_en_vec)
        new_perm = bz_grid.new_sort_perm(1, 3, 0, w_en, w_vc)
        # each row of new_perm should be range(n_br) if the sorting is stable
        test_stability = np.array([(x == range(n_br)).all() for x in new_perm])
        self.assertTrue(test_stability.all())
        test_result = np.isclose(en_vec, s_en_vec)
        if not test_result.all():
            np.set_printoptions(linewidth=2000)
            np.set_printoptions(threshold=sys.maxsize)
            for (t_r, a_v, b_v) in zip(test_result, en_vec, s_en_vec):
                if (~t_r).any():
                    t_p = (~t_r).any(0)
                    print(np.concatenate((a_v[:, t_p], b_v[:, t_p]), axis=1))
            print("\nThere are ",
                  (~test_result).any(2).any(1).sum(),
                  " grid points with swapped values\n")
        self.assertTrue(test_result.all())


if __name__ == '__main__':
    unittest.main()
