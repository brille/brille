#!/usr/bin/env python3
import os
import sys
import unittest
import numpy as np
from pathlib import Path

addpaths = [Path(), Path('..'), Path('../..')] # need outer nest too
config = os.environ.get('CMAKE_CONFIG_TYPE') # set by ctest -C <cfg>
if config:
  for path in addpaths:
    if Path(path, config).exists():
      addpaths.append(Path(path,config))
sys.path[:0] = [str(path.absolute()) for path in addpaths]

import brille 

def get_latvec(lens, angs):
    # Calculates the lattice vectors after the convention of self.spglib
    xhat = np.array([1, 0, 0]);
    yhat = np.array([np.cos(angs[2]), np.sin(angs[2]), 0]);
    vol = np.sqrt(1 - np.sum(np.cos(angs)**2) + 2*np.prod(np.cos(angs)));
    zhat = np.array([np.cos(angs[1]), \
                     (np.cos(angs[0]) - np.cos(angs[1])*np.cos(angs[2]))/np.sin(angs[2]), \
                     vol / np.sin(angs[2])]);
    return np.array([lens[0]*xhat, lens[1]*yhat, lens[2]*zhat]);


class UtilsTestBZ (unittest.TestCase):
    # Tests create_bz routines
    a = 4.0
    b = 4.0
    c = 5.0
    alpha = np.pi / 2
    beta = np.pi / 2
    gamma = 2 * np.pi / 3
    spg = 'P 6'

    def check_lattice(self, bz, vals=None):
        lattice = bz.lattice.star
        if not vals:
            vals = [self.a, self.b, self.c, self.alpha, self.beta, self.gamma]
        self.assertAlmostEqual(lattice.a, vals[0])
        self.assertAlmostEqual(lattice.b, vals[1])
        self.assertAlmostEqual(lattice.c, vals[2])
        self.assertAlmostEqual(lattice.alpha, vals[3])
        self.assertAlmostEqual(lattice.beta, vals[4])
        self.assertAlmostEqual(lattice.gamma, vals[5])

    def test_separate_latpars(self):
        bz = brille.utils.create_bz(self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.spg)
        self.check_lattice(bz)
        # Again with keywords
        bz = brille.utils.create_bz(a=self.a, b=self.b, c=self.c, alpha=self.alpha, beta=self.beta, gamma=self.gamma, spacegroup=self.spg)
        self.check_lattice(bz)

    def test_separate_latpars_degrees(self):
        alf, bet, gam = np.array([self.alpha, self.beta, self.gamma]) / np.pi * 180
        bz = brille.utils.create_bz(self.a, self.b, self.c, alf, bet, gam, self.spg)
        self.check_lattice(bz)

    def test_lengths_angles(self):
        bz = brille.utils.create_bz([self.a, self.b, self.c], [self.alpha, self.beta, self.gamma], self.spg)
        self.check_lattice(bz)
        # With keywords
        bz = brille.utils.create_bz(lens=[self.a, self.b, self.c], angs=[self.alpha, self.beta, self.gamma], spacegroup=self.spg)
        self.check_lattice(bz)

    def test_lengths_angles_degrees(self):
        bz = brille.utils.create_bz([self.a, self.b, self.c], np.array([self.alpha, self.beta, self.gamma]) / np.pi * 180, self.spg)
        self.check_lattice(bz)

    def test_lattice_vectors(self):
        latvec = get_latvec([self.a, self.b, self.c], [self.alpha, self.beta, self.gamma])
        bz = brille.utils.create_bz(latvec, self.spg)
        self.check_lattice(bz)

    def test_invalid_input_wrong_spg(self):
        # Tests incorrect spacegroup for lattice parameters - cannot find IR wedge
        with self.assertRaises(RuntimeError):
            bz = brille.utils.create_bz([4, 5, 6], [90, 90, 90], 'P 6')
        bz = brille.utils.create_bz([4, 5, 6], [90, 90, 90], 'P 6', wedge_search=False)
        self.check_lattice(bz, [4, 5, 6] + [np.pi/2]*3)
        
    def test_invalid_input_negative_vol(self):
        # Tests for physically impossible lattice
        with self.assertRaises(ValueError):
            bz = brille.utils.create_bz([4, 5, 6], [30, 30, 90], 'P 6')

    def test_no_spacegroup(self):
        # Tests if no spacegroup is given
        with self.assertRaises(ValueError):
            bz = brille.utils.create_bz([4, 5, 6], [90, 90, 90])

    def test_spacegroup_number(self):
        # Tests for a spacegroup / Hall number instead of a Hall symbol
        bz = brille.utils.create_bz(self.a, self.b, self.c, self.alpha, self.beta, self.gamma, 1)
        self.check_lattice(bz)

    def test_invalid_spacegroups(self):
        # Tests for invalid spacegroups given
        with self.assertRaises(ValueError):
            # Must be a string or number
            bz = brille.utils.create_bz([4, 5, 6], [90, 90, 90], [])
        with self.assertRaises(RuntimeError):
            # "Unknown lattice type"
            bz = brille.utils.create_bz([4, 5, 6], [90, 90, 90], -1)

    def test_using_reciprocal(self):
        # Test with reciprocal lattice
        latvec = get_latvec([self.a, self.b, self.c], [self.alpha, self.beta, self.gamma])
        bz = brille.utils.create_bz(2 * np.pi * np.linalg.inv(latvec).T, self.spg, is_reciprocal=True)
        self.check_lattice(bz)


class UtilsTestGrid (unittest.TestCase):
    # Tests create_grid routines
    bz = brille.utils.create_bz(4, 4, 6, 90, 90, 110, 'P 2')

    def check_type(self, expected_type, *args, **kwargs):
        obj = brille.utils.create_grid(self.bz, *args, **kwargs)
        self.assertTrue(isinstance(obj, expected_type))

    def test_expected_types(self):
        # Checks the input arguments results in the correct grid types
        self.check_type(brille.BZTrellisQdd, node_volume_fraction=0.1)
        self.check_type(brille.BZTrellisQdc, node_volume_fraction=0.1, complex_vectors=True)
        self.check_type(brille.BZTrellisQcc, node_volume_fraction=0.1, complex_vectors=True, complex_values=True)
        self.check_type(brille.BZMeshQdd, mesh=True)
        self.check_type(brille.BZMeshQdc, mesh=True, complex_vectors=True)
        self.check_type(brille.BZMeshQcc, mesh=True, complex_vectors=True, complex_values=True)
        self.check_type(brille.BZNestQdd, nest=True, max_volume=0.1)
        self.check_type(brille.BZNestQdc, nest=True, max_volume=0.1, complex_vectors=True)
        self.check_type(brille.BZNestQcc, nest=True, max_volume=0.1, complex_vectors=True, complex_values=True)
        self.check_type(brille.BZNestQdd, nest=True, number_density=10.0)
        self.check_type(brille.BZNestQdc, nest=True, number_density=10.0, complex_vectors=True)
        self.check_type(brille.BZNestQcc, nest=True, number_density=10.0, complex_vectors=True, complex_values=True)

    def test_bad_input(self):
        # Checks for invalid input
        with self.assertRaises(ValueError):
            # Specifying both mesh and nest [cannot decide]
            brille.utils.create_grid(self.bz, mesh=True, nest=True)
        with self.assertRaises(ValueError):
            # Parameters for BZNestQ* without nest=True
            brille.utils.create_grid(self.bz, max_volume=0.1)
        with self.assertRaises(ValueError):
            # Missing arguments `max_volume` or `number_density`
            brille.utils.create_grid(self.bz, nest=True) 


if __name__ == '__main__':
  unittest.main()
