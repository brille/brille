#!/usr/bin/env python3
import unittest
import numpy as np
from pathlib import Path

from load_local import load
br_mod = load(('_brille', 'brille._brille'), prefer_installed=True, search=[Path(), Path('..')])
br_py = load(('brille.utils',), prefer_installed=True, search=[Path(), Path('..'), Path('../..')])

if br_mod != br_py:
    print(f"Loaded binary module {br_mod.__file__} and python module {br_py.__file__} differ")

def get_latvec(lens, angs):
    # Calculates the lattice vectors after the convention of self.spglib
    xhat = np.array([1, 0, 0])
    yhat = np.array([np.cos(angs[2]), np.sin(angs[2]), 0])
    vol = np.sqrt(1 - np.sum(np.cos(angs)**2) + 2*np.prod(np.cos(angs)))
    zhat = np.array([np.cos(angs[1]),
                     (np.cos(angs[0]) - np.cos(angs[1])*np.cos(angs[2]))/np.sin(angs[2]),
                     vol / np.sin(angs[2])])
    return np.array([lens[0]*xhat, lens[1]*yhat, lens[2]*zhat])


class UtilsTestBZ (unittest.TestCase):
    # Tests create_bz routines
    a, b, c, alpha, beta, gamma = 4.0, 4.0, 5.0, np.pi/2, np.pi/2, 2*np.pi/3
    spg = 'P 6'

    def check_lattice(self, bz, vals=None):
        lattice = bz.lattice
        if not vals:
            vals = [self.a, self.b, self.c, self.alpha, self.beta, self.gamma]
        self.assertAlmostEqual(lattice.a, vals[0])
        self.assertAlmostEqual(lattice.b, vals[1])
        self.assertAlmostEqual(lattice.c, vals[2])
        self.assertAlmostEqual(lattice.alpha, vals[3])
        self.assertAlmostEqual(lattice.beta, vals[4])
        self.assertAlmostEqual(lattice.gamma, vals[5])

    def test_separate_latpars(self):
        bz = br_py.utils.create_bz(self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.spg)
        self.check_lattice(bz)
        # Again with keywords
        bz = br_py.utils.create_bz(a=self.a, b=self.b, c=self.c, alpha=self.alpha, beta=self.beta, gamma=self.gamma, spacegroup=self.spg)
        self.check_lattice(bz)

    def test_separate_latpars_degrees(self):
        alf, bet, gam = np.array([self.alpha, self.beta, self.gamma]) / np.pi * 180
        bz = br_py.utils.create_bz(self.a, self.b, self.c, alf, bet, gam, self.spg)
        self.check_lattice(bz)

    def test_lengths_angles(self):
        bz = br_py.utils.create_bz([self.a, self.b, self.c], [self.alpha, self.beta, self.gamma], self.spg)
        self.check_lattice(bz)
        # With keywords
        bz = br_py.utils.create_bz(lens=[self.a, self.b, self.c], angs=[self.alpha, self.beta, self.gamma], spacegroup=self.spg)
        self.check_lattice(bz)

    def test_lengths_angles_degrees(self):
        bz = br_py.utils.create_bz([self.a, self.b, self.c], np.array([self.alpha, self.beta, self.gamma]) / np.pi * 180, self.spg)
        self.check_lattice(bz)

    def test_lattice_vectors(self):
        latvec = get_latvec([self.a, self.b, self.c], [self.alpha, self.beta, self.gamma])
        bz = br_py.utils.create_bz(latvec, self.spg)
        self.check_lattice(bz)

    def test_invalid_input_wrong_spg(self):
        # Tests incorrect spacegroup for lattice parameters - cannot find IR wedge
        with self.assertRaises(RuntimeError):
            bz = br_py.utils.create_bz([4, 5, 6], [90, 90, 90], 'P 6')
        bz = br_py.utils.create_bz([4, 5, 6], [90, 90, 90], 'P 6', wedge_search=False)
        self.check_lattice(bz, [4, 5, 6] + [np.pi/2]*3)
        
    def test_invalid_input_negative_squared_vol(self):
        # Tests for physically impossible lattice
        with self.assertRaises(ValueError):
            bz = br_py.utils.create_bz([4, 5, 6], [30, 30, 90], 'P 6')

    def test_no_spacegroup(self):
        # Tests if no spacegroup is given
        with self.assertRaises(ValueError):
            bz = br_py.utils.create_bz([4, 5, 6], [90, 90, 90])

    # 'Hall number' specification removed in v0.6.0
    # def test_spacegroup_number(self):
    #     # Tests for a spacegroup / Hall number instead of a Hall symbol
    #     bz = br_py.utils.create_bz(self.a, self.b, self.c, self.alpha, self.beta, self.gamma, 1)
    #     self.check_lattice(bz)

    def test_invalid_spacegroups(self):
        # Tests for invalid spacegroups given
        with self.assertRaises(ValueError):
            # Must be a string
            bz = br_py.utils.create_bz([4, 5, 6], [90, 90, 90], [])
        with self.assertRaises(ValueError):
            # Must be a string
            bz = br_py.utils.create_bz([4, 5, 6], [90, 90, 90], -1)

    def test_using_reciprocal(self):
        # Test with reciprocal lattice
        latvec = get_latvec([self.a, self.b, self.c], [self.alpha, self.beta, self.gamma])
        bz = br_py.utils.create_bz(2 * np.pi * np.linalg.inv(latvec).T, self.spg, is_reciprocal=True)
        self.check_lattice(bz)


class UtilsTestGrid (unittest.TestCase):
    # Tests create_grid routines
    bz = br_py.utils.create_bz(4, 4, 6, 90, 90, 110, 'P 2')

    def check_type(self, expected_type, *args, **kwargs):
        obj = br_py.utils.create_grid(self.bz, *args, **kwargs)
        self.assertTrue(isinstance(obj, expected_type))

    def test_expected_types(self):
        # Checks the input arguments results in the correct grid types
        self.check_type(br_mod.BZTrellisQdd, node_volume_fraction=0.1)
        self.check_type(br_mod.BZTrellisQdc, node_volume_fraction=0.1, complex_vectors=True)
        self.check_type(br_mod.BZTrellisQcc, node_volume_fraction=0.1, complex_vectors=True, complex_values=True)
        self.check_type(br_mod.BZMeshQdd, mesh=True)
        self.check_type(br_mod.BZMeshQdc, mesh=True, complex_vectors=True)
        self.check_type(br_mod.BZMeshQcc, mesh=True, complex_vectors=True, complex_values=True)
        self.check_type(br_mod.BZNestQdd, nest=True, max_volume=0.1)
        self.check_type(br_mod.BZNestQdc, nest=True, max_volume=0.1, complex_vectors=True)
        self.check_type(br_mod.BZNestQcc, nest=True, max_volume=0.1, complex_vectors=True, complex_values=True)
        self.check_type(br_mod.BZNestQdd, nest=True, number_density=10.0)
        self.check_type(br_mod.BZNestQdc, nest=True, number_density=10.0, complex_vectors=True)
        self.check_type(br_mod.BZNestQcc, nest=True, number_density=10.0, complex_vectors=True, complex_values=True)

    def test_bad_input(self):
        # Checks for invalid input
        with self.assertRaises(ValueError):
            # Specifying both mesh and nest [cannot decide]
            br_py.utils.create_grid(self.bz, mesh=True, nest=True)
        with self.assertRaises(ValueError):
            # Parameters for BZNestQ* without nest=True
            br_py.utils.create_grid(self.bz, max_volume=0.1)
        with self.assertRaises(ValueError):
            # Missing arguments `max_volume` or `number_density`
            br_py.utils.create_grid(self.bz, nest=True)


if __name__ == '__main__':
  unittest.main()
