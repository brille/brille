"""Normalization of interpolated eigenvectors (set_vector_normalization)."""
import numpy as np
import pytest
from brille import BrillouinZone, BZMeshQdc, Lattice

GAMMA = 2  # RotatesLike: phonon eigenvectors
ANGSTROM, REAL_LATTICE = 1, 3  # LengthUnit values


def grid():
    # one atom, so brille can transform the vectors as phonon eigenvectors (Gamma)
    lattice = Lattice(((3, 3, 3), (90, 90, 90)), spacegroup="P 1", basis=([[0, 0, 0]], [0]))
    bz = BrillouinZone(lattice)
    return BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 200)


def unit_vectors(q):
    """Two branches of unit vectors that turn quickly with q, so interpolation shortens them."""
    a = 2 * np.pi * q[:, 0]
    phase = np.exp(2j * np.pi * q[:, 1])
    b0 = np.stack([np.cos(a), np.sin(a), np.zeros_like(a)], axis=1) * phase[:, None]
    b1 = np.stack([np.zeros_like(a), np.cos(3 * a), np.sin(3 * a)], axis=1)
    return np.ascontiguousarray(np.stack([b0, b1], axis=1))  # (n, 2 branches, 3)


def filled(g=None):
    g = grid() if g is None else g
    q = np.asarray(g.rlu)
    # phonon-like eigenvectors in Cartesian units, as Euphonic provides them
    g.fill(np.ones((len(q), 2)), (1,), unit_vectors(q), (0, 3, 0, GAMMA, ANGSTROM))
    return g


def queries():
    return np.random.default_rng(3).random((200, 3)) - 0.5


def norms(g, metric=None):
    _, vecs = g.ir_interpolate_at(queries())
    v = np.asarray(vecs).reshape(len(queries()), 2, 3)
    m = np.ones(3) if metric is None else np.asarray(metric)
    return np.einsum("qbi,i,qbi->qb", v.conj(), m, v).real


def test_off_by_default_and_interpolation_shortens_vectors():
    g = filled()
    assert not g.normalizes_vectors
    n = norms(g)
    assert n.max() <= 1 + 1e-12
    assert n.min() < 0.99  # the interpolated vectors really are short


def test_normalization_gives_unit_vectors():
    g = filled()
    g.set_vector_normalization()
    assert g.normalizes_vectors
    assert list(g.vector_metric) == []
    np.testing.assert_allclose(norms(g), 1, atol=1e-12)


def test_metric_normalization_keeps_the_sign():
    g = filled()
    eta = [1.0, 1.0, -1.0]
    before = norms(g, eta)
    g.set_vector_normalization(metric=eta)
    after = norms(g, eta)
    np.testing.assert_allclose(np.abs(after), 1, atol=1e-12)
    assert np.all(np.sign(after) == np.sign(before))


def test_metric_must_match_the_branch_length():
    g = filled()
    with pytest.raises(RuntimeError, match="one weight for each of the 3 elements"):
        g.set_vector_normalization(metric=[1.0, 1.0])


def test_turning_normalization_off_again():
    g = filled()
    g.set_vector_normalization()
    g.set_vector_normalization(False)
    assert norms(g).min() < 0.99


def test_setting_survives_fill():
    g = grid()
    g.set_vector_normalization()
    filled(g)
    np.testing.assert_allclose(norms(g), 1, atol=1e-12)


def test_setting_survives_hdf5(tmp_path):
    g = filled()
    g.set_vector_normalization(metric=[1.0, 1.0, -1.0])
    path = str(tmp_path / "grid.h5")
    g.to_file(path)
    loaded = BZMeshQdc.from_file(path)
    assert loaded.normalizes_vectors
    assert list(loaded.vector_metric) == [1.0, 1.0, -1.0]
    np.testing.assert_allclose(norms(loaded, [1, 1, -1]), norms(g, [1, 1, -1]), atol=1e-12)


def test_lattice_units_are_refused():
    """A vector's length in lattice units depends on the lattice, which the interpolator doesn't know."""
    g = grid()
    q = np.asarray(g.rlu)
    g.fill(np.ones((len(q), 2)), (1,), unit_vectors(q), (0, 3, 0, GAMMA, REAL_LATTICE))
    with pytest.raises(RuntimeError, match="Cartesian units"):
        g.set_vector_normalization()
    # and refilling a normalizing grid with lattice-unit vectors fails before changing it
    g = filled()
    g.set_vector_normalization()
    with pytest.raises(RuntimeError, match="Cartesian units"):
        g.fill(np.ones((len(q), 2)), (1,), unit_vectors(q), (0, 3, 0, GAMMA, REAL_LATTICE))
    np.testing.assert_allclose(norms(g), 1, atol=1e-12)

