#!/usr/bin/env python3
import unittest
import numpy as np

def fetchLoad(loader, fetchfile, **kwds):
    """Fetch remote binary file from the :py:mod:`brille` repository

    Parameters
    ----------
    loader : function
        The function which loads the requested file into Python
    fetchfile : str
        The filename (with extension) of the file to fecth
    **kwds :
        All keyword arguments are passed to the `loader` function

    Returns
    -------
    loaded
        The loaded object or value returned by the `loader` function
    """
    import tempfile
    import requests
    import pathlib
    base_url = "https://raw.githubusercontent.com/brille/brille/master/wrap/tests"
    with tempfile.TemporaryDirectory() as tmp_dir:
        r = requests.get(base_url + "/" + fetchfile)
        if not r.ok:
            raise Exception("Fetching {} failed with reason '{}'".format(fetchfile, r.reason))
        out_path = str(pathlib.Path(tmp_dir, fetchfile))
        open(out_path, 'wb').write(r.content)
        return loader(out_path, **kwds)


def getLoad(loader, file, **kwds):
    """Load a binary file

    If the file is not located in the `brille` module directory this function
    will attempt to obtain it from the `brille` repository over any available
    network connection.

    Parameters
    ----------
    loader : function
        The function which loads the requested file into Python
    fetchfile : str
        The filename (with extension) of the file to fecth
    **kwds :
        All keyword arguments are passed to the `loader` function

    Returns
    -------
    loaded
        The loaded object or value returned by the `loader` function
    """
    import os
    import pathlib
    docs_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        return loader(str(pathlib.Path(docs_dir, file)), **kwds)
    except FileNotFoundError:
        print('{} not found in {}. Fetching remote content.'.format(file, docs_dir))
        return fetchLoad(loader, file, **kwds)


class GammaTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Set up nacl test data, lattice and grid

        from brille import Basis, Symmetry, Lattice, BrillouinZone, BZTrellisQdc
        # load numpy arrays from the compressed binary pack
        nacl = getLoad(np.load, 'test_5_gamma.npz')
        # construct the NaCl direct lattice
        bas = Basis(nacl['atom_positions'], nacl['atom_index'])
        sym = Symmetry(nacl['spacegroup_mat'], nacl['spacegroup_vec'])
        lat = Lattice(nacl['basis_vectors'], symmetry=sym, basis=bas)

        # use it to construct an irreducible Brillouin zone
        bz = BrillouinZone(lat)
        # and use that to produce a hybrid interpolation grid
        # with parameters stored in the binary pack
        max_volume = float(nacl['grid_max_volume'])
        always_triangulate = bool(nacl['grid_always_triangulate'])
        grid = BZTrellisQdc(bz, max_volume, always_triangulate)

        cls.nacl = nacl
        cls.lat = lat
        cls.grid = grid

    def test_nacl(self):
        # Verify interpolation with defined test data/lattice/grid
        nacl = self.nacl
        lat = self.lat
        grid = self.grid

        # verify the stored grid points to ensure we have the same irreducible
        # wedge and grid
        # Though the points could be permuted:
        perm = np.hstack([np.argwhere(np.all(np.isclose(nacl['grid_rlu'], x), axis=1)) for x in grid.rlu]).flatten()
        self.assertTrue(np.allclose(nacl['grid_rlu'][perm], grid.rlu))

        # insert the Euphonic-derived eigenvalues and eigenvectors for the grid
        # points, which are used in the interpolation
        # Original vectors elements enumeration was [0, 24 0, 3, 0, 0] where
        # RotatesLike::Gamma = 3. New interface includes extra enumeration
        # LengthUnit::real_lattice = 3, and the RotatesLike enumeration has
        # changed, so RotatesLike::Gamma = 2 now. Replace here rather than
        # regenerate file.
        vec_els = np.array([0, 24, 0, 2, 3, 0, 0])
        grid.fill(
            nacl['grid_values'][perm], nacl['grid_values_elements'], nacl['grid_values_weights'],
            nacl['grid_vectors'][perm], vec_els, nacl['grid_vectors_weights'],
            bool(nacl['grid_sort']))

        # the fourth grid point is inside of the irreducible volume, which
        # ensures we avoid degeneracies
        q_ir = nacl['grid_rlu'][4:5]
        # find all symmetry equivalent q within the first Brillouin zone by
        # applying each pointgroup operator to q_ir to find q_nu = R^T_nu q_ir
        q_nu = np.einsum('xji,aj->xi', lat.pointgroup.W, q_ir)

        # The std::sort algorithm does not provide the same pointgroup sorting
        # on all systems, but there should be a permutation mapping:
        perm = np.array([np.squeeze(np.argwhere(np.all(np.isclose(q_nu, x), axis=1))) for x in nacl['q_nu']])
        # The permutation must be complete and unique
        self.assertTrue(np.unique(perm).size == q_nu.shape[0])
        # And the permuted q_nu must match the stored q_nu
        q_nu = q_nu[perm]
        self.assertTrue(np.allclose(q_nu, nacl['q_nu']))

        # Use the grid to interpolate at each q_nu:
        br_val, br_vec = grid.ir_interpolate_at(q_nu)
        br_val = np.squeeze(br_val)

        # verify that q_ir does not have degeneracies:
        self.assertFalse(np.any(np.isclose(np.diff(br_val, axis=1), 0.)))
        # and that the 'interpolated' eigenvalues are identical for all q_nu
        self.assertTrue(np.allclose(np.diff(br_val, axis=0), 0.))
        # plus that the interpolated eigenvalues match the store Euphonic eigenvalues
        self.assertTrue(np.allclose(br_val, nacl['euphonic_values']))

        # convert the eigenvectors into the same cartesian coordinate system
        # used by Euphonic
        br_vec = np.einsum('ba,ijkb->ijka', nacl['basis_vectors'], br_vec)
        # load the Euphonic calculated eigenvectors
        eu_vec = nacl['euphonic_vectors']
        # The 'interpolated' eigenvectors and the Euphonic eigenvectors should
        # only be equivalent up to an overall phase factor, so find it:
        antiphase = np.exp(-1J * np.angle(np.einsum('qmij,qmij->qm', np.conj(eu_vec), br_vec)))
        # and remove the phase from the interpolated eigenvectors
        br_vec = np.einsum('ab,abij->abij', antiphase, br_vec)
        # now all eigenvectors must match
        self.assertTrue(np.allclose(br_vec, eu_vec))

        # as an extra check we *could* verfiy that the output of brille is
        # unchanged, but we don't care about changes in overall phase factor

    def test_cartesian_eigenvector_interpolation(self):
        # Verify interpolation with Cartesian eigenvectors
        nacl = self.nacl
        grid = self.grid

        # Permuted grid points:
        perm = np.hstack([
            np.argwhere(np.all(np.isclose(nacl['grid_rlu'], x), axis=1))
                for x in grid.rlu]
            ).flatten()

        # Convert input grid vectors from basis to Cartesian
        grid_vecs_cart =  np.einsum(
            'ba,ijkb->ijka', nacl['basis_vectors'], nacl['grid_vectors'])
        # Use RotatesLike::Gamma = 2 and LengthUnit::angstrom = 1
        vec_els = np.array([0, 24, 0, 2, 1, 0, 0])
        grid.fill(
            nacl['grid_values'][perm], nacl['grid_values_elements'],
            nacl['grid_values_weights'], grid_vecs_cart[perm],
            vec_els, nacl['grid_vectors_weights'])

        # Use the grid to interpolate at each q_nu:
        br_val, br_vec = grid.ir_interpolate_at(nacl['q_nu'])

        # load the Euphonic calculated eigenvectors
        eu_vec = nacl['euphonic_vectors']
        # The 'interpolated' eigenvectors and the Euphonic eigenvectors should
        # only be equivalent up to an overall phase factor, so find it:
        antiphase = np.exp(-1J * np.angle(
            np.einsum('qmij,qmij->qm', np.conj(eu_vec), br_vec)))
        # and remove the phase from the interpolated eigenvectors
        br_vec = np.einsum('ab,abij->abij', antiphase, br_vec)
        # now all eigenvectors must match
        self.assertTrue(np.allclose(br_vec, eu_vec))

    def test_unsupported_rotateslike_lengthunit_combinations_error(self):
        nacl = self.nacl
        grid = self.grid
        # RotatesLike, Lengthunit
        combos = [(0, 1), (1, 0), (2, 2), (2, 4)]

        for combo in combos:
            vec_els = np.array([0, 24, 0, *combo, 0, 0])
            grid.fill(nacl['grid_values'], nacl['grid_values_elements'],
            nacl['grid_values_weights'], nacl['grid_vectors'],
            vec_els, nacl['grid_vectors_weights'])
            with self.assertRaises(RuntimeError):
                grid.ir_interpolate_at(nacl['q_nu'])


if __name__ == '__main__':
    unittest.main()

#
# # The data used in  the above test was generated using NDLT1145, a Dell XPS
# # laptop running Windows 10, using
# #   brille   v0.5.0-master.1b5e5a6
# #   brilleu  v0.2.3
# #   euphonic v0.3.0
# # and the following script:
#
# import brille as b, numpy as np, euphonic, brilleu
# from brilleu.brilleu import getBrillEuObj
#
# nacl = {}
# nacl['brille_version'] = b.version
# nacl['brilleu_version'] = brilleu.__version__
# nacl['euphonic_version'] = euphonic.__version__
# nacl['grid_sort'] = False
# nacl['grid_max_volume'] = 0.1
#
# breu = getBrillEuObj('NaCl', sort=nacl['grid_sort'], max_volume=nacl['grid_max_volume'], parallel=False)
#
# nacl['basis_vectors'] = breu.grid.BrillouinZone.lattice.star.lattice_matrix
# nacl['atom_positions'] = breu.crystal.atom_positions
# nacl['atom_index'] = breu.crystal.atom_index
# nacl['spacegroup_mat'] = breu.crystal.symmetry.W
# nacl['spacegroup_vec'] = breu.crystal.symmetry.w
#
# nacl['grid_always_triangulate'] = False
# nacl['grid_rlu'] = breu.grid.rlu
# nacl['grid_values'] = breu.grid.values
# nacl['grid_values_elements'] = (1.,)
# nacl['grid_values_weights'] = (13605.693, 0., 0.)
# nacl['grid_vectors'] = breu.grid.vectors
# nacl['grid_vectors_elements'] = (0, 3*breu.data.crystal.n_atoms, 0, 3, 0, 0)
# nacl['grid_vectors_weights'] = (0., 1., 0.)
#
# q_ir = breu.grid.rlu[4:5]
# nacl['pointgroup'] = breu.grid.BrillouinZone.lattice.pointgroup.W
# nacl['q_nu'] = np.einsum('xji,aj->xi', nacl['pointgroup'], q_ir)
#
# br_nu = breu.QpointPhononModes(nacl['q_nu'], interpolate=True)
# nacl['brille_values'] = br_nu.ω
# nacl['brille_vectors'] = br_nu.ε
#
# eu_nu = breu.QpointPhononModes(nacl['q_nu'], interpolate=False)
# nacl['euphonic_values'] = eu_nu.ω
# nacl['euphonic_vectors'] = eu_nu.ε
#
# np.savez('test_5_gamma.npz', **nacl)
