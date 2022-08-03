#!/usr/bin/env python3
import unittest
import numpy as np
import brille as br_py
from brille import _brille as br_mod


class Lattice(unittest.TestCase):
    def setUp(self) -> None:
        self.cmo_vectors = np.array([[1.7797736800000001, 1.027552813333333, 6.219738366666666],
                                     [-1.7797736800000001, 1.027552813333333, 6.219738366666666],
                                     [0.0, -2.055105626666667, 6.219738366666666]])
        self.cmo_positions = np.array([[0.89401258, 0.89401259, 0.89401258],
                                       [0.10598742, 0.10598741, 0.10598742],
                                       [0.5, 0.5, 0.5],
                                       [0., 0.99999999, 0.]])
        self.cmo_types = [0, 0, 1, 2]
        self.cmo_rotations = np.array([[[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
                                       [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
                                       [[0, 0, -1], [0, -1, 0], [-1, 0, 0]]])
        self.cmo_translations = np.array([[0., 0., 0.] for _ in self.cmo_rotations])

        self.cmo_avg_positions = np.array([[0.8940125833333333, 0.8940125833333333, 0.8940125833333333],
                                           [0.1059874166666667, 0.1059874166666666, 0.1059874166666667],
                                           [0.5, 0.5, 0.5],
                                           [0., 1., 0.]])

    def test_a_init2(self):
        a = 1.0
        b = np.pi / 2
        v = 1.0
        l = br_mod.Lattice([a, a, a], [b, b, b])  # compute the volume when required
        self.assertEqual(l.a, a)
        self.assertEqual(l.b, a)
        self.assertEqual(l.c, a)
        self.assertEqual(l.alpha, b)
        self.assertEqual(l.beta, b)
        self.assertEqual(l.gamma, b)
        self.assertAlmostEqual(l.volume, v)  # since 0.9999999999999999 != 1.0

    def test_a_init3(self):
        a = (1, 1, 1)
        b = np.array([1, 1, 1]) * np.pi / 2
        v = 1.0
        l = br_mod.Lattice(a, b)
        self.assertEqual(l.a, a[0])
        self.assertEqual(l.b, a[1])
        self.assertEqual(l.c, a[2])
        self.assertEqual(l.alpha, b[0])
        self.assertEqual(l.beta, b[1])
        self.assertEqual(l.gamma, b[2])
        self.assertAlmostEqual(l.volume, v)  # since 0.9999999999999999 != 1.0

    def test_b_tensors(self):
        a = (6, 6, 9)
        b = np.array([90, 90, 120]) * np.pi / 180
        l = br_mod.Lattice(a, b)
        cmt = np.array([[a[0] ** 2, a[0] * a[1] * np.cos(b[2]), 0],
                        [a[0] * a[1] * np.cos(b[2]), a[1] ** 2, 0],
                        [0, 0, a[2] ** 2]])
        self.assertAlmostEqual((l.get_covariant_metric_tensor() - cmt).sum(), 0)
        cnt = np.linalg.inv(cmt)
        self.assertAlmostEqual((l.get_contravariant_metric_tensor() - cnt).sum(), 0)

    def test_d_equality(self):
        a = 1.0
        b = np.pi / 2
        l1 = br_mod.Lattice([a, a, a], [b, b, b])
        l2 = br_mod.Lattice([a, a, a], [b, b, b])
        self.assertEqual(l1, l2)

    def test_e_duals(self):
        d = br_mod.Lattice((1, 1, 1), np.array([1, 1, 1]) * np.pi / 2, real_space=True)
        r = br_mod.Lattice(np.array([1, 1, 1]) * np.pi * 2, np.array([1, 1, 1]) * np.pi / 2, real_space=False)
        self.assertEqual(d, r)

    def test_f_basis(self):
        basis = br_mod.Basis(self.cmo_positions, self.cmo_types)
        generators = br_mod.Symmetry(self.cmo_rotations, self.cmo_translations)
        symmetry = generators.generate()

        # Without specifying snap_to_symmetry the basis positions are as specified
        lat = br_mod.Lattice(self.cmo_vectors, symmetry, basis, snap_to_symmetry=False)
        self.assertTrue(np.allclose(lat.basis.positions, self.cmo_positions))

        # but specifying snap_to_symmetry=True forces the positions to take on their symmetry averages
        lat = br_mod.Lattice(self.cmo_vectors, symmetry, basis, snap_to_symmetry=True)
        self.assertTrue(np.allclose(lat.basis.positions, self.cmo_avg_positions))

    def test_g_wrapped_constructor(self):
        # build tuples here just to avoid typing the same thing multiple times
        bt = self.cmo_positions, self.cmo_types
        st = self.cmo_rotations, self.cmo_translations

        # Providing a Basis with no symmetry information should raise a ValueError
        self.assertRaises(ValueError, br_py.Lattice, self.cmo_vectors, basis=bt)

        # Providing both a spacegroup string (or strings) and Symmetry information should raise an error
        self.assertRaises(ValueError, br_py.Lattice, self.cmo_vectors, spacegroup='ITName', symmetry=st)
        self.assertRaises(ValueError, br_py.Lattice, self.cmo_vectors, spacegroup=('HMSymbol', 'HMChoice'), symmetry=st)

        # Providing Symmetry and Basis as positional arguments raises an error
        self.assertRaises(TypeError, br_py.Lattice, self.cmo_vectors, br_mod.Symmetry(*st), br_mod.Basis(*bt))

        # Until we can stop using Python 3.6 and 3.7 the first argument *can not* be position-only
        # # The single required argument can not be provided as a keyword
        # self.assertRaises(TypeError, br_py.Lattice, values=self.cmo_vectors)

        # And none of kwargs can have the same name 'values'
        self.assertRaises(TypeError, br_py.Lattice, self.cmo_vectors, values=True)

        # Providing Symmetry and Basis as keyword arguments works
        lat = br_py.Lattice(self.cmo_vectors, symmetry=br_mod.Symmetry(*st), basis=br_mod.Basis(*bt))
        self.assertTrue(np.allclose(lat.basis.positions, self.cmo_positions))

        # The objects can be created for us from tuple arguments if we wish:
        lat = br_py.Lattice(self.cmo_vectors, symmetry=st, basis=bt)
        self.assertTrue(np.allclose(lat.basis.positions, self.cmo_positions))

        # and of course, we can pass extra keywords (and their order doesn't matter)
        lat = br_py.Lattice(self.cmo_vectors, basis=bt, snap_to_symmetry=False, symmetry=st)
        self.assertTrue(np.allclose(lat.basis.positions, self.cmo_avg_positions))


if __name__ == '__main__':
    unittest.main()
