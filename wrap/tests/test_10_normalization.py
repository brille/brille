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


def filled_in_lattice_units():
    g = grid()
    q = np.asarray(g.rlu)
    g.fill(np.ones((len(q), 2)), (1,), unit_vectors(q), (0, 3, 0, GAMMA, REAL_LATTICE))
    return g


def test_cartesian_eigenvectors_are_normalized_by_default():
    g = filled()
    assert g.vector_normalization == "automatic"
    assert g.normalizes_vectors
    assert list(g.vector_metric) == []
    np.testing.assert_allclose(norms(g), 1, atol=1e-12)


def test_without_normalization_interpolation_shortens_vectors():
    g = filled()
    g.set_vector_normalization(False)
    assert g.vector_normalization == "off"
    assert not g.normalizes_vectors
    n = norms(g)
    assert n.max() <= 1 + 1e-12
    assert n.min() < 0.99  # the interpolated vectors really are short


def test_lattice_units_are_not_normalized_by_default():
    """A vector's length in lattice units depends on the lattice, which the interpolator doesn't know."""
    g = filled_in_lattice_units()
    assert g.vector_normalization == "automatic"
    assert not g.normalizes_vectors
    assert norms(g).min() < 0.99  # left as interpolated, and no error


def test_forcing_normalization_refuses_lattice_units():
    g = filled_in_lattice_units()
    with pytest.raises(RuntimeError, match="Cartesian units"):
        g.set_vector_normalization(True)
    # and refilling a force-normalized grid with lattice-unit vectors fails before changing it
    g = filled()
    g.set_vector_normalization(True)
    q = np.asarray(g.rlu)
    with pytest.raises(RuntimeError, match="Cartesian units"):
        g.fill(np.ones((len(q), 2)), (1,), unit_vectors(q), (0, 3, 0, GAMMA, REAL_LATTICE))
    np.testing.assert_allclose(norms(g), 1, atol=1e-12)


def test_none_restores_the_automatic_default():
    g = filled()
    g.set_vector_normalization(False)
    g.set_vector_normalization(None)
    assert g.vector_normalization == "automatic"
    np.testing.assert_allclose(norms(g), 1, atol=1e-12)


def test_metric_normalization_keeps_the_sign():
    g = filled()
    eta = [1.0, 1.0, -1.0]
    g.set_vector_normalization(False)
    before = norms(g, eta)
    g.set_vector_normalization(metric=eta)
    after = norms(g, eta)
    np.testing.assert_allclose(np.abs(after), 1, atol=1e-12)
    assert np.all(np.sign(after) == np.sign(before))


def test_metric_must_match_the_branch_length():
    g = filled()
    with pytest.raises(RuntimeError, match="one weight for each of the 3 elements"):
        g.set_vector_normalization(metric=[1.0, 1.0])


def test_setting_survives_fill():
    g = grid()
    g.set_vector_normalization(False)
    filled(g)
    assert g.vector_normalization == "off"
    assert norms(g).min() < 0.99


@pytest.mark.parametrize("normalize, metric, mode", [
    (None, None, "automatic"), (True, [1.0, 1.0, -1.0], "on"), (False, None, "off")])
def test_setting_survives_hdf5(tmp_path, normalize, metric, mode):
    g = filled()
    g.set_vector_normalization(normalize, metric)
    path = str(tmp_path / "grid.h5")
    g.to_file(path)
    loaded = BZMeshQdc.from_file(path)
    assert loaded.vector_normalization == mode
    assert list(loaded.vector_metric) == (metric or [])
    np.testing.assert_allclose(norms(loaded), norms(g), atol=1e-12)
