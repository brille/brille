"""Moving points into the first Brillouin zone and the irreducible wedge.

The per-point search runs on plain numbers and in parallel; these check its
invariants on several lattices, including centred ones (moved in the primitive
cell) and time-reversal symmetry, and that thread count does not change results.
"""
import numpy as np
import pytest
from brille import BrillouinZone, Lattice

CASES = {
    "fcc": (((5.64, 5.64, 5.64), (90, 90, 90)), "-F 4 2 3", False),
    "hexagonal": (((3.0, 3.0, 5.0), (90, 90, 120)), "-P 6 2", False),
    "rhombohedral": (((4.9, 4.9, 13.8), (90, 90, 120)), "-R 3 2\"", False),
    "monoclinic C": (((5.1, 8.8, 5.2), (90, 104.5, 90)), "-C 2y", False),
    "triclinic": (((4.02, 4.90, 3.29), (98.97, 86.24, 88.47)), "-P 1", False),
    "P3, time reversal": (((4.1, 4.1, 5.3), (90, 90, 120)), "P 3", True),
}


def points(bz, n=4000):
    rng = np.random.default_rng(5)
    return np.concatenate([(rng.random((n, 3)) - 0.5) * 8,
                           rng.integers(-12, 13, (n, 3)) / rng.choice([2, 3, 4, 6], (n, 1)),
                           np.asarray(bz.vertices)])


@pytest.mark.parametrize("name", CASES)
def test_moveinto(name):
    lp, hall, tr = CASES[name]
    bz = BrillouinZone(Lattice(lp, spacegroup=hall), time_reversal_symmetry=tr)
    Q = points(bz)
    q, tau = bz.moveinto(Q, threads=1)
    q, tau = np.asarray(q), np.asarray(tau)
    np.testing.assert_array_equal(tau, np.round(tau))
    np.testing.assert_allclose(q + tau, Q, atol=1e-12)
    assert all(bz.isinside(q))
    q4, tau4 = bz.moveinto(Q, threads=4)
    np.testing.assert_array_equal(np.asarray(q4), q)
    np.testing.assert_array_equal(np.asarray(tau4), tau)


def inside_ir_polyhedron(bz, q, tol=1e-9):
    """Whether each q is in the irreducible polyhedron.

    Containment in a convex polyhedron survives a linear map, so the faces'
    planes can be found in lattice units with plain vector arithmetic.
    """
    vertices = np.asarray(bz.ir_vertices)
    centre = vertices.mean(axis=0)
    inside = np.ones(len(q), dtype=bool)
    for face in bz.ir_vertices_per_face:
        polygon = vertices[face]
        # Newell's normal is robust to nearly collinear vertices
        normal = np.cross(polygon, np.roll(polygon, -1, axis=0)).sum(axis=0)
        normal /= np.linalg.norm(normal)
        point = polygon.mean(axis=0)
        if np.dot(normal, centre - point) > 0:
            normal = -normal
        inside &= (q - point) @ normal <= tol * np.abs(vertices).max()
    return inside


@pytest.mark.parametrize("name", CASES)
def test_ir_moveinto(name):
    lp, hall, tr = CASES[name]
    bz = BrillouinZone(Lattice(lp, spacegroup=hall), time_reversal_symmetry=tr)
    Q = points(bz)
    q, tau, R, invR = (np.asarray(x) for x in bz.ir_moveinto(Q, threads=1))
    np.testing.assert_allclose(np.einsum("nji,nj->ni", R, q) + tau, Q, atol=1e-10)
    np.testing.assert_array_equal(np.einsum("nij,njk->nik", R, invR), np.broadcast_to(np.eye(3, dtype=int), R.shape))
    assert all(bz.isinside(q))
    assert all(inside_ir_polyhedron(bz, q))
    q4, tau4, R4, invR4 = (np.asarray(x) for x in bz.ir_moveinto(Q, threads=4))
    for a, b in ((q, q4), (tau, tau4), (R, R4), (invR, invR4)):
        np.testing.assert_array_equal(a, b)
