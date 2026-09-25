"""Lattices given as basis vectors with explicit symmetry operations.

Phonon codes give a primitive cell as vectors, from spglib or similar, and
the vectors carry round-off: an off-diagonal 4e-16 in a tetragonal cell, say.
brille used to take that metric as it was, and search for the irreducible wedge
with tolerances. Some centred lattices then got no irreducible zone, and some
primitive ones got a zone whose mesh ignored ``max_size``. Now the metric is
made invariant under the point group, and the wedge is the point group's
Dirichlet cone, whose planes have integer coefficients.

``explicit_lattices.json`` holds, by ITA number, a primitive cell from the
reference-model harness (AFLOW lattice, spglib reduction) with its operations.
"""
import json
from pathlib import Path

import numpy as np
import pytest
from brille import BrillouinZone, BZMeshQdc, Lattice

CASES = json.loads((Path(__file__).parent / "explicit_lattices.json").read_text())


def lattice(number):
    case = CASES[number]
    return Lattice(np.array(case["lattice"]), symmetry=(np.array(case["rotations"]), np.array(case["translations"])))


@pytest.mark.parametrize("number", CASES)
def test_irreducible_zone_and_mesh(number):
    bz = BrillouinZone(lattice(number))
    operations = len({tuple(np.ravel(r)) for r in CASES[number]["rotations"]})
    assert bz.polyhedron.volume / bz.ir_polyhedron.volume == pytest.approx(operations)
    coarse, fine = (BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / k).rlu.shape[0] for k in (100, 3000))
    assert fine > 5 * coarse > 5 * 4, "the mesh must refine as max_size shrinks"


@pytest.mark.parametrize("number", CASES)
def test_metric_is_invariant(number):
    lat = lattice(number)
    metric = np.array(lat.get_covariant_metric_tensor())
    for r in np.array(CASES[number]["rotations"]):
        assert np.allclose(r.T @ metric @ r, metric, rtol=0, atol=1e-14 * np.abs(metric).max())


def test_vectors_keep_their_orientation():
    """Only round-off is removed; the Cartesian frame of the vectors stays"""
    given = np.array(CASES["127"]["lattice"])
    assert given[1, 0] != 0  # the round-off this test is about
    lat = lattice("127")
    vectors = np.array(lat.real_vectors).T  # given as rows, returned as columns
    assert np.allclose(vectors, given, rtol=0, atol=1e-12)
    metric = np.array(lat.get_covariant_metric_tensor())
    assert metric[0, 1] == metric[0, 2] == metric[1, 2] == 0
    assert metric[0, 0] == metric[1, 1]


def test_rotated_vectors_keep_their_orientation():
    """A cell in an arbitrary Cartesian frame stays in that frame"""
    angle = 0.3
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    rotation = rotation @ np.array([[1, 0, 0], [0, np.cos(1.1), -np.sin(1.1)], [0, np.sin(1.1), np.cos(1.1)]])
    case = CASES["221"]
    given = np.array(case["lattice"]) @ rotation.T
    lat = Lattice(given, symmetry=(np.array(case["rotations"]), np.array(case["translations"])))
    assert np.allclose(np.array(lat.real_vectors).T, given, rtol=0, atol=1e-12)
    metric = np.array(lat.get_covariant_metric_tensor())
    assert metric[0, 0] == metric[1, 1] == metric[2, 2]
    assert metric[0, 1] == metric[0, 2] == metric[1, 2] == 0
