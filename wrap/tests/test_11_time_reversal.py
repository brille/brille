"""Time-reversal symmetry for phonon eigenvectors in crystals without inversion.

Time reversal is anti-unitary: e(-q) = e(q)*. Combined with a real operation R,
the eigenvector at -Rq is the complex conjugate of the one at Rq. brille used to
add time reversal as a space inversion, which has no atom mapping in a crystal
without inversion, and refused.

Both Rq and -Rq are reconstructed from the same interpolated e(q) in the
irreducible zone, so the relation holds exactly for any stored data.
"""
import numpy as np
import pytest
from brille import BrillouinZone, BZMeshQdc, Lattice

GAMMA, ANGSTROM = 2, 1  # RotatesLike, LengthUnit


def filled_grid(spacegroup, positions, types):
    lattice = Lattice(((4.0, 5.0, 6.0), (90, 100, 90)), spacegroup=spacegroup, basis=(positions, types))
    bz = BrillouinZone(lattice, time_reversal_symmetry=True)
    grid = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 100)
    n_atoms, n_modes = len(positions), 3 * len(positions)
    rng = np.random.default_rng(7)
    n = len(grid.rlu)
    vectors = rng.standard_normal((n, n_modes, 3 * n_atoms)) + 1j * rng.standard_normal((n, n_modes, 3 * n_atoms))
    values = rng.random((n, n_modes))
    grid.fill(values, (1,), np.ascontiguousarray(vectors), (0, 3 * n_atoms, 0, GAMMA, ANGSTROM))
    return bz, grid


def interior_point(bz):
    return np.asarray(bz.ir_polyhedron.vertices).mean(axis=0)


@pytest.mark.parametrize("spacegroup, positions, types, rotation", [
    # P1: time reversal is paired with the identity
    ("P 1", [[0.0, 0.0, 0.0], [0.3, 0.2, 0.1]], [0, 1], np.eye(3)),
    # P2 (unique axis b): the 2-fold swaps the two atoms, so its atom mapping matters
    ("P 2y", [[0.1, 0.2, 0.3], [-0.1, 0.2, -0.3]], [0, 0], np.diag([-1, 1, -1])),
])
def test_minus_rq_is_the_conjugate_of_rq(spacegroup, positions, types, rotation):
    bz, grid = filled_grid(spacegroup, positions, types)
    q = interior_point(bz)
    rq = rotation @ q
    values, vectors = grid.ir_interpolate_at(np.array([rq, -rq]))
    values, vectors = np.asarray(values), np.asarray(vectors)
    np.testing.assert_allclose(values[1], values[0], atol=1e-12)
    np.testing.assert_allclose(vectors[1], np.conj(vectors[0]), atol=1e-12)
