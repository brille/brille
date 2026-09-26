"""The C++ grid of the structured mesh against the Python reference.

A lattice's Delaunay triangulation doesn't depend on implementation choices: it is
unique, or, where a Selling parameter is zero, its cells are split by a symmetric
rule. So the C++ grid must equal the reference tetrahedron for tetrahedron,
exactly (grid points are small rationals).
"""
import numpy as np
import pytest
from brille import Lattice
from brille._brille import _LatticeGrid
from reference.lattice_tri.grid import Grid as ReferenceGrid, D

CASES = {
    "triclinic": (((4.02, 4.90, 3.29), (98.97, 86.24, 88.47)), "-P 1"),
    "monoclinic P": (((5.1, 8.8, 5.2), (90, 104.5, 90)), "-P 2y"),
    "orthorhombic P": (((3.1, 4.4, 5.3), (90, 90, 90)), "-P 2 2"),
    "tetragonal P": (((4.0, 4.0, 6.1), (90, 90, 90)), "-P 4 2"),
    "hexagonal": (((3.0, 3.0, 5.0), (90, 90, 120)), "-P 6 2"),
    "cubic P": (((4.0,) * 3, (90,) * 3), "-P 4 2 3"),
}
# primitive cells, whose reciprocal lattices are not degenerate (symmetry P 1 only)
PRIMITIVE = {
    "fcc primitive": np.array([[0, 2.82, 2.82], [2.82, 0, 2.82], [2.82, 2.82, 0]]),
    "rhombohedral primitive": np.array([[2.45, 1.41451, 4.6], [-2.45, 1.41451, 4.6], [0, -2.82902, 4.6]]),
}


def grid_input(name):
    """(reciprocal basis rows, point group on reciprocal lattice coordinates)"""
    if name in PRIMITIVE:
        lattice = Lattice(PRIMITIVE[name], spacegroup="P 1")
    else:
        lattice = Lattice(CASES[name][0], spacegroup=CASES[name][1])
    basis = np.asarray(lattice.reciprocal_vectors, dtype=float).T          # rows
    ops = [np.asarray(w, dtype=np.int64).T for w in np.asarray(lattice.pointgroup.W)]
    return basis, ops


NAMES = list(CASES) + list(PRIMITIVE)


@pytest.mark.parametrize("name", NAMES)
def test_grid_matches_reference(name):
    basis, ops = grid_input(name)
    grid, reference = _LatticeGrid(basis, ops), ReferenceGrid(basis, ops)
    assert grid.degenerate == reference.degenerate
    radius = 1.2 * np.linalg.norm(basis, axis=1).max()
    factor = D // _LatticeGrid.scale
    def inner(tets, points_of):
        out = set()
        for t in tets:
            pts = points_of(t)
            if np.linalg.norm(np.mean([np.asarray(p, float) / D @ basis for p in pts], axis=0)) < 0.6 * radius:
                out.add(frozenset(pts))
        return out
    ours = inner(grid.patch(radius), lambda t: [tuple(int(c) * factor for c in p) for p in t])
    theirs = inner(reference.patch(radius), lambda t: list(t))
    assert ours and ours == theirs


@pytest.mark.parametrize("name", NAMES)
def test_grid_is_invariant(name):
    basis, ops = grid_input(name)
    assert _LatticeGrid(basis, ops).invariant(2 * np.linalg.norm(basis, axis=1).max())


@pytest.mark.parametrize("name", ["triclinic", "hexagonal", "fcc primitive"])
def test_grid_locates_points(name):
    basis, ops = grid_input(name)
    grid = _LatticeGrid(basis, ops)
    rng = np.random.default_rng(3)
    for x in rng.uniform(-3, 3, (500, 3)):
        tet, weights = grid.locate(x)
        corners = np.asarray(tet[0], float) / _LatticeGrid.scale @ basis
        assert min(weights) > -1e-12
        np.testing.assert_allclose(np.asarray(weights) @ corners, x, atol=1e-10)
