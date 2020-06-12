#!/usr/bin/env python3
"""Validate the phonon eigenvector rotation method for correctness"""
import os
import sys
import unittest
import tempfile
import requests
try:
    from importlib.util import find_spec
except ImportError:
    raise Exception("This test requires importlib.util (and Python3)")
import numpy as np
# Now the actual search for the module
if find_spec('brille') is None or find_spec('brille._brille') is None:
    raise Exception('Required compiled brille module not found')
if find_spec('euphonic') is None or find_spec('euphonic.data') is None or find_spec('euphonic.data.interpolation') is None:
    raise Exception("Required euphonic.data.interpolation module not found!")
import brille as b
from pathlib import Path
from euphonic.data.interpolation import InterpolationData
from brille.euphonic import BrEu

def fetch_InterpolationData_object(material):
    base_url = "https://raw.githubusercontent.com/g5t/brille/master/docs/tutorials"
    file_to_fetch = material + ".castep_bin"
    with tempfile.TemporaryDirectory() as tmp_dir:
        r = requests.get(base_url + "/" + file_to_fetch)
        if not r.ok:
            raise Exception("Fetching {} failed with reason '{}'".format(file_to_fetch, r.reason))
        out_path = Path(tmp_dir, file_to_fetch)
        open(str(out_path), 'wb').write(r.content)
        idata = InterpolationData.from_castep(path=tmp_dir, seedname=material)
    return idata

def get_InterpolationData_object(material):
	validation_dir = os.path.dirname(os.path.abspath(__file__))
    tutorial_dir = str(Path(validation_dir, '..', 'docs', 'tutorials'))
	try:
		return InterpolationData.from_castep(path=tutorial_dir, seedname=material)
	except FileNotFoundError:
		print('{} not found in {}. Fetching remote content.'.format(material, tutorial_dir))
	return fetch_InterpolationData_object(material)

class PhononEigenvectorRotationValidator(unittest.TestCase):
    """A TestCase object class to validate phonon eigenvector rotation by brille"""

    def test_NaCl(self):
        # fetch and load the NaCl.castep_bin file from the brille repository
        idata = get_InterpolationData_object('NaCl')
        # Do not sort the modes on neighbouring trellis vertices to make comparison with Euphonic easier
        breu = BrEu(idata, sort=False, trellis=True, max_volume=0.1, parallel=False)
        # pick a trellis Q vertex inside the irreducible polyhedron, away from the boundaries:
        q_ir = breu.grid.rlu[6:7] # shape = (1,3)
        # pull together the pointgroup operations
        ptgr = b.PointSymmetry(breu.grid.BrillouinZone.lattice.hall)
        # and apply each to q_ir to find q_nu = R^T_nu q_ir
        q_nu = np.einsum('xji,aj->xi', ptgr.W, q_ir)

        # use brille to find eigenenergies and eigenvectors at each q_nu
        # this *should* be only an application of the rotation, so
        #   all omega_nu are expected to be identical
        #   and all epsilon_nu will be permuted as Gamma(q|nu) dictates
        br_omega_nu, br_epsilon_nu = breu.frqs_vecs(q_nu, interpolate=True)
        # use Euphonic to diagonalise the dynamical matrix at each q_nu
        eu_omega_nu, eu_epsilon_nu = breu.frqs_vecs(q_nu, interpolate=False)

        # verify that all brille eigenvalues are identical for each q_nu
        self.assertTrue(np.allclose(np.diff(br_omega_nu.magnitude, axis=0), 0.))
        # and that these results match the Euphonic results
        self.assertTrue(np.allclose(br_omega_nu.magnitude, eu_omega_nu.magnitude))
        # Now that we're sure the permutation is the same, we can verify the eigenvectors
        # brille stores and returns eigenvectors expressed in units of the conventional direct lattice
        # while euphonic calculates, returns, and stores them in a cartesian coordinate system defined
        # relative to the lattice by its lat_vec matrix.
        # A brille.spglib method handles conversion between these descriptions:
        br_epsilon_nu_xyz = breu.brspgl.conventional_to_orthogonal_eigenvectors(br_epsilon_nu)

        # The eigenvectors returned by brille and those from Euphonic should be the same up to an overall complex phase factor
        # the eigenvectors are normalised, so the phase is calculated directly
        antiphase_per_q_per_mode = np.exp(-1J*np.angle(np.einsum('qmij,qmij->qm', np.conj(eu_epsilon_nu), br_epsilon_nu_xyz)))
        # and removed from the cartesian coordinate eigenvectors
        br_epsilon_nu_xyz_phased = np.array([[x0*y0 for x0,y0 in zip(x,y)] for x,y in zip(antiphase_per_q_per_mode, br_epsilon_nu_xyz)])
        # now all eigenvectors returned from brille should be the same as
        # calculated by Euphonic directly
        self.assertTrue(np.allclose(br_epsilon_nu_xyz_phased, eu_epsilon_nu))
