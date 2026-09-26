"""The structured mesh (the lattice grid clipped to the irreducible zone) has the
properties any mesh of the irreducible zone must have (mesh_properties.py): it
fills the zone exactly, it is conforming, and its boundary matches itself under
the zone's face pairings, which the TetGen mesh does not (GitHub #114, test_18).
"""
import itertools

import numpy as np
import pytest
from brille import Lattice
from brille._brille import _LatticeTri
import mesh_properties as mp
from test_20_lattice_boundary import CASES, CENTRING, primitive_input


def primitive_basis(name):
    lengths_angles, symmetry, centring, _ = CASES[name]
    A = np.asarray(Lattice(lengths_angles, spacegroup=symmetry).real_vectors, dtype=float).T
    return 2 * np.pi * np.linalg.inv(CENTRING[centring] @ A).T


def properties(name, n):
    metric, ops = primitive_input(name)
    B = primitive_basis(name)
    mesh = _LatticeTri(metric, ops, n)
    vertices = np.asarray(mesh.vertices) @ B
    tetrahedra = np.asarray(mesh.tetrahedra)
    faces = [np.asarray(f) @ B for f in mesh.faces]
    centre = np.concatenate(faces).mean(axis=0)
    planes = []
    for polygon in faces:
        normal = np.linalg.svd(polygon - polygon.mean(axis=0))[2][-1]
        d = normal @ polygon.mean(axis=0)
        if normal @ centre > d:
            normal, d = -normal, -d
        planes.append((normal, d))
    operations = [B.T @ np.asarray(o, float) @ np.linalg.inv(B.T) for o in ops]
    translations = [np.array(s) @ B for s in itertools.product(range(-2, 3), repeat=3)]
    tol = 1e-8 * np.abs(vertices).max()
    return dict(volume=mp.total_volume(vertices, tetrahedra) / (abs(np.linalg.det(B)) / len(ops)),
                open_faces=mp.open_faces(vertices, tetrahedra, planes, tol),
                boundary=mp.boundary_mismatches(vertices, tetrahedra, planes, operations, translations, tol))


@pytest.mark.parametrize("n", (2, 3))
@pytest.mark.parametrize("name", CASES)
def test_structured_mesh_properties(name, n):
    result = properties(name, n)
    assert result["volume"] == pytest.approx(1, abs=1e-9)
    assert result["open_faces"] == 0
    checked, missing = result["boundary"]
    assert checked > 0 and missing == 0, f"{missing} of {checked} boundary images missing"
