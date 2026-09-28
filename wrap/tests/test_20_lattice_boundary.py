"""The C++ boundary of the irreducible zone against the Python reference.

The irreducible polyhedron, its pairing cells and its special points are fixed by
geometry, not by implementation choices (e.g. tie-breaks), so the C++ and the
exact-arithmetic reference must agree: the same faces, the same cells on each
face, and the same special points, to a tolerance for the C++'s coordinates.
"""
import numpy as np
import pytest
from brille import Lattice
from brille._brille import _LatticeBoundary
from reference.lattice_tri.boundary import boundary as reference_boundary

CENTRING = {
    "P": np.eye(3),
    "C": np.array([[0.5, 0.5, 0], [-0.5, 0.5, 0], [0, 0, 1]]),
    "I": np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
    "F": np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
    "R": np.array([[2 / 3, 1 / 3, 1 / 3], [-1 / 3, 1 / 3, 1 / 3], [-1 / 3, -2 / 3, 1 / 3]]),
}
CASES = {
    "triclinic": (((4.02, 4.90, 3.29), (98.97, 86.24, 88.47)), "-P 1", "P", True),
    "monoclinic C": (((5.1, 8.8, 5.2), (90, 104.5, 90)), "-C 2y", "C", False),
    "tetragonal I": (((4.0, 4.0, 9.0), (90, 90, 90)), "-I 4 2", "I", False),
    "trigonal P 3, time reversal": (((3.1, 3.1, 4.7), (90, 90, 120)), "P 3", "P", True),
    "hexagonal 6/m": (((3.0, 3.0, 5.0), (90, 90, 120)), "-P 6", "P", False),
    "rhombohedral -3m": (((4.9, 4.9, 13.8), (90, 90, 120)), '-R 3 2"', "R", False),
    "fcc": (((5.64,) * 3, (90,) * 3), "-F 4 2 3", "F", False),
}


def primitive_input(name):
    """(reciprocal metric, point group) in the primitive reciprocal basis"""
    lengths_angles, symmetry, centring, time_reversal = CASES[name]
    lattice = Lattice(lengths_angles, spacegroup=symmetry)
    A = np.asarray(lattice.real_vectors, dtype=float).T              # rows: conventional direct
    P = CENTRING[centring] @ A                                        # rows: primitive direct
    B = 2 * np.pi * np.linalg.inv(P).T                                # rows: primitive reciprocal
    W = np.asarray(lattice.pointgroup.W)
    cart = [A.T @ w @ np.linalg.inv(A.T) for w in W]                  # Cartesian rotations
    ops = []
    for r in cart + ([-r for r in cart] if time_reversal else []):
        m = np.linalg.inv(B.T) @ r @ B.T
        mi = np.round(m).astype(np.int64)
        assert np.allclose(m, mi, atol=1e-8)
        if not any((mi == o).all() for o in ops):
            ops.append(mi)
    return B @ B.T, ops


def as_set(points, digits=9):
    return frozenset(tuple(np.round(np.asarray(p, float), digits) + 0.0) for p in points)


@pytest.mark.parametrize("name", CASES)
def test_boundary_matches_reference(name):
    metric, ops = primitive_input(name)
    ours = _LatticeBoundary(metric, ops)
    planes, face_vertices, cells, special = reference_boundary(metric, ops)
    theirs_faces = {as_set([[float(c) for c in v] for v in fv]): i for i, fv in enumerate(face_vertices)}
    ours_faces = {as_set(f): i for i, f in enumerate(ours.faces)}
    assert set(ours_faces) == set(theirs_faces)
    for key, i in ours_faces.items():
        j = theirs_faces[key]
        ours_cells = sorted(sorted(as_set(c)) for c in ours.cells[i])
        theirs_cells = sorted(sorted(as_set([[float(x) for x in v] for v in c])) for c in cells[j])
        assert ours_cells == theirs_cells
    assert as_set(ours.special_points) == as_set([[float(x) for x in p] for p in special])
