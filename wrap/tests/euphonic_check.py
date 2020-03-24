#!/usr/bin/env python3
# This file should not be renamed 'test*.py' until module/submodule loading is
# figured out from within a build directory.
# At the moment, if you've built the _brille binary but not installed it for
# python to find on its path, this script puts the empty (__import__.py only)
# brille 'package' on the python search path in order to load the euphonic
# submodule. When the submodule is loaded the empty 'package' gets loaded as
# well, and then the whole thing fails to run.
"""Run tests of BrEu and Euphonic."""
import os
import sys
import unittest
from euphonic.data.interpolation import InterpolationData
try:
    from importlib.util import find_spec
except ImportError:
    raise Exception("This test requires importlib.util (and Python3)")
import numpy as np

# Try to find the brille module:
# It might be in the current directory, or a sub-directory called, e.g., Debug
ADDPATH = os.getcwd()
if os.path.exists('Debug'):
    ADDPATH = os.path.join(ADDPATH, 'Debug')
sys.path.append(ADDPATH)
# Make sure we can load the binary module before the euphonic submodule tries:
if find_spec('brille') is not None and find_spec('brille._brille') is not None:
    import brille
elif find_spec('_brille') is not None:
    # pylint: disable=e0401
    import _brille as brille
else:
    raise Exception("Required brille module not found!")
# We need to find the pure-python submodule brille.euphonic:
sys.path.append(os.path.split(os.getcwd())[0])
if find_spec('brille') is not None and find_spec('brille.euphonic') is not None:
    from brille.euphonic import BrEu
else:
    raise Exception("Required brille.euphonic module not found!")


def load_interpolation_data(named):
    """Load a data file from the tests folder."""
    b = find_spec('brille')
    if b is None:
        bdir = '.' # punt! maybe we're in a build directory?
    else:
        bdir = b.submodule_search_locations[0]
    testdir = os.path.join(bdir, '..', 'tests')
    if not os.path.exists(testdir):
        raise Exception('Could not locate the tests directory')
    return InterpolationData(seedname=named, path=testdir)


def hermitian_product(v_0, v_1, first=None, last=None):
    """Find the hermitian product between two sets of vectors."""
    return np.dot(np.conj(v_0[first:last].flatten()), v_1[first:last].flatten())


class TestEuphonic(unittest.TestCase):
    """A TestCase object class to run tests of the BrEu object."""

    def test_a(self):
        """Test whether SimPhony gives the same result *at* grid points.

        This test should always pass, since the BZGridQcomplex object is filled
        by calculating eigen-energies and eigen-vectors using SimPhony.
        Unfortunately, due to rounding errors, this test can fail to run if an
        out-of-bounds interpolation point is required by brille.

        The moveinto keyword is what makes it possible to verify the
        "interpolation" result against the SimPhony result. This same keyword
        also opens up the possibility of rounding-errors leading to a runtime
        error of brille.
        With the moveinto keyword omitted (or set to True), out-of-zone grid
        points *will* likely have different polarization vectors compared to
        the SimPhony result as one interpolates at q=Q-τ and the other
        calculates at Q. It's unclear if this difference is important for
        intensity calculations.
        """
        i_data = load_interpolation_data('nb')
        symsim = BrEu(i_data, halfN=(2, 2, 2))

        q_rlu = symsim.grid.rlu
        int_freq, int_vecs = symsim.frqs_vecs(q_rlu, interpolate=True,
                                              moveinto=False)
        sim_freq, sim_vecs = symsim.frqs_vecs(q_rlu, interpolate=False)

        ad_freq = np.abs(int_freq - sim_freq).magnitude
        as_freq = np.abs(int_freq + sim_freq).magnitude
        # Check if |interpolated-simulated|/|interpolated+simulated|>1e-14
        # AND if |interpolated-simulated|>1e-14
        unequal_freq = (ad_freq > 1e-14*as_freq) * (ad_freq > 1e-14)
        self.assertFalse(unequal_freq.any())

        # The vectors only need to be equal up to an arbitrary phase
        # which is equivalent to saying that the inner product between
        # equal ion eigenvectors for each branch should have ||ϵ⋅ϵ||²≡1
        n_pt, n_br, n_io, n_d = int_vecs.shape
        int_vecs = int_vecs.reshape(n_pt*n_br*n_io, n_d)
        sim_vecs = sim_vecs.reshape(n_pt*n_br*n_io, n_d)
        product = [hermitian_product(x, y) for x, y in zip(int_vecs, sim_vecs)]

        self.assertTrue(np.isclose(np.abs(product), 1).all())

    def test_b(self):
        """Do something."""
        i_data = load_interpolation_data('nb')
        symsim = BrEu(i_data, halfN=(10, 10, 10))
        print(symsim.grid.centre_sort_perm())


if __name__ == '__main__':
    unittest.main()
