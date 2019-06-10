#!/usr/bin/env python3
"""Run tests of SymSim and Simphony."""
import os
import sys
import unittest
from simphony.data.interpolation import InterpolationData
try:
    from importlib.util import find_spec
except ImportError:
    raise Exception("This test requires importlib.util (and Python3)")
import numpy as np

# Try to find the symbz module:
# It might be in the current directory, or a sub-directory called, e.g., Debug
ADDPATH = os.getcwd()
if os.path.exists('Debug'):
    ADDPATH = os.path.join(ADDPATH, 'Debug')
sys.path.append(ADDPATH)
# Now the actual search for the module
if find_spec('symbz') is not None and find_spec('symbz._symbz') is not None:
    import symbz
elif find_spec('_symbz') is not None:
    import _symbz as symbz
else:
    raise Exception("Required symbz module not found!")
# We need to find the pure-python submodule symbz.simphony:
sys.path.append(os.path.split(os.getcwd())[0])
if find_spec('symbz') is not None and find_spec('symbz.simphony') is not None:
    from symbz.simphony import SymSim
else:
    raise Exception("Required symbz.simphony module not found!")

def load_interpolation_data(named):
    test_spec = find_spec('tests')
    if test_spec is None:
        raise Exception('Could not locate the tests module directory')
    seed = os.path.join(test_spec.submodule_search_locations[0], named)
    return InterpolationData(seed)

class TestSimphony(unittest.TestCase):
    """A TestCase object class to run tests of the SymSim object."""

    def test_a(self):
        """Test whether SimPhony gives the same result *at* grid points.

        This test should always pass, since the BZGridQcomplex object is filled
        by calculating eigen-energies and eigen-vectors using SimPhony.
        Unfortunately, due to rounding errors, this test can fail to run if an
        out-of-bounds interpolation point is required by symbz.

        The moveinto keyword is what makes it possible to verify the
        "interpolation" result against the SimPhony result. This same keyword
        also opens up the possibility of rounding-errors leading to a runtime
        error of symbz.
        With the moveinto keyword omitted (or set to True), out-of-zone grid
        points *will* likely have different polarization vectors compared to
        the SimPhony result as one interpolates at q=Q-τ and the other
        calculates at Q. It's unclear if this difference is important for
        intensity calculations.
        """
        i_data = load_interpolation_data('nb')
        symsim = SymSim(i_data, halfN=(2, 2, 2))

        q_rlu = symsim.grid.mapped_rlu
        int_freq, int_vecs = symsim.frqs_vecs(q_rlu, interpolate=True, moveinto=False)
        sim_freq, sim_vecs = symsim.frqs_vecs(q_rlu, interpolate=False)

        ad_freq = np.abs(int_freq - sim_freq).magnitude
        as_freq = np.abs(int_freq + sim_freq).magnitude
        # Check if |interpolated-simulated|/|interpolated+simulated|>1e-14
        # AND if |interpolated-simulated|>1e-14
        unequal_freq = (ad_freq > 1e-14*as_freq) * (ad_freq > 1e-14)
        self.assertFalse(unequal_freq.any())

        ad_vecs = np.abs(int_vecs - sim_vecs)
        as_vecs = np.abs(int_vecs + sim_vecs)
        unequal_vecs = (ad_vecs > 1e-14*as_vecs) * (ad_vecs > 1e-14)
        self.assertFalse(unequal_vecs.any())

if __name__ == '__main__':
    unittest.main()
