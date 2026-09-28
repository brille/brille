"""conventional_to_primitive recovers a primitive cell's modes from its conventional cell's.

A spring model defined in real space is solved on both cells: the conventional
cell's modes at q, converted, must be the primitive cell's modes at q, for the
atoms brille keeps (Lattice.primitive_basis), in the cell phase convention.
"""
import itertools

import numpy as np
import pytest
from brille import Lattice
from brille.utils import conventional_to_primitive

SIX_P = {"I": [-3, 3, 3, 3, -3, 3, 3, 3, -3], "F": [0, 3, 3, 3, 0, 3, 3, 3, 0],
         "C": [3, 3, 0, -3, 3, 0, 0, 0, 6], "R": [4, -2, -2, 2, 2, -4, 2, 2, 2]}
CASES = {
    "F": ("Fm-3m", ([5.6] * 3, [90] * 3), [[0, 0, 0], [0.5, 0.5, 0.5]]),
    "I": ("Im-3m", ([3.0] * 3, [90] * 3), [[0, 0, 0]]),
    "C": ("C 1 2/m 1", ([5.0, 6.0, 4.0], [90, 104, 90]), [[0, 0, 0], [0.3, 0, 0.2], [0.7, 0, 0.8]]),
    "R": ("R-3m", ([3.8, 3.8, 9.5], [90, 90, 120]), [[0, 0, 0], [0, 0, 0.5]]),
}


def cartesian(lengths, angles):
    a, b, c = lengths
    al, be, ga = np.radians(angles)
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    return np.array([[a, 0, 0], [b * np.cos(ga), b * np.sin(ga), 0], [cx, cy, np.sqrt(c * c - cx * cx - cy * cy)]])


def dynamical_matrix(q_frac, rows, positions, types):
    """Cell-convention dynamical matrix of central springs, q in the cell's reciprocal lattice units"""
    n, cutoff = len(positions), 7.0
    reach = cutoff + 2 * np.abs(positions @ rows).max() + 1
    extent = [int(np.ceil(reach * np.linalg.norm(v))) + 1 for v in np.linalg.inv(rows).T]
    cells = np.array(list(itertools.product(*(range(-e, e + 1) for e in extent))))
    D = np.zeros((3 * n, 3 * n), dtype=complex)
    for i in range(n):
        for j in range(n):
            d = (positions[j] - positions[i] + cells) @ rows
            r = np.linalg.norm(d, axis=1)
            near = (r > 1e-9) & (r < cutoff)
            k = (1.0 + 0.4 * (types[i] + types[j])) * np.exp(-(r[near] / 3.0) ** 2)
            phi = -np.einsum("p,pa,pb->pab", k / r[near] ** 2, d[near], d[near])
            D[3 * i:3 * i + 3, 3 * j:3 * j + 3] += np.einsum("pab,p->ab", phi, np.exp(2j * np.pi * cells[near] @ q_frac))
            D[3 * i:3 * i + 3, 3 * i:3 * i + 3] -= phi.sum(axis=0)
    mass = np.repeat(1.0 + np.asarray(types, dtype=float), 3)
    return (D / np.sqrt(np.outer(mass, mass)) + (D / np.sqrt(np.outer(mass, mass))).conj().T) / 2


def solve(q_frac, rows, positions, types):
    w2, e = np.linalg.eigh(dynamical_matrix(q_frac, rows, positions, types))
    return w2, e.T.reshape(len(w2), len(positions), 3)


@pytest.mark.parametrize("centring", CASES)
def test_conversion_matches_the_primitive_cell(centring):
    symbol, cell, sites = CASES[centring]
    lattice = Lattice(cell, spacegroup=symbol)
    shifts = np.asarray(lattice.centring_vectors)
    positions = np.vstack([(np.asarray(sites) + t) % 1 for t in shifts])
    types = np.tile(np.minimum(np.arange(len(sites)), 1), len(shifts))
    lattice = Lattice(cell, spacegroup=symbol, basis=(positions, types))
    # the positions as brille holds them, which it may have fitted to the symmetry
    positions, types = np.asarray(lattice.basis.positions), np.asarray(lattice.basis.types)
    kept = np.asarray(lattice.primitive_basis.positions)
    kept_types = np.asarray(lattice.primitive_basis.types)
    rows = cartesian(*cell)
    P = np.array(SIX_P[centring], dtype=float).reshape(3, 3) / 6
    q = np.array([[0, 0, 0], [0.13, 0.27, 0.41], [0.5, 0, 0], [1, 0, 0], [0.25, 0.25, 0.25]])
    conventional = [solve(x, rows, positions, types) for x in q]
    values, vectors = conventional_to_primitive(
        lattice, q, np.array([c[0] for c in conventional]), np.array([c[1] for c in conventional]))
    for i, x in enumerate(q):
        # the primitive cell: its vectors are P's columns, its atoms those brille kept
        w2, e = solve(x @ P, P.T @ rows, kept @ np.linalg.inv(P).T, kept_types)
        assert np.allclose(np.sort(values[i]), w2, atol=1e-9)
        for value in np.unique(np.round(w2, 8)):
            mine = vectors[i][np.abs(values[i] - value) < 1e-7].reshape(-1, 3 * len(kept))
            exact = e[np.abs(w2 - value) < 1e-7].reshape(-1, 3 * len(kept))
            assert np.abs(mine.T @ mine.conj() - exact.T @ exact.conj()).max() < 1e-8
