"""The structured mesh (the lattice grid clipped to the irreducible zone) has the
properties any mesh of the irreducible zone must have (mesh_properties.py): it
fills the zone exactly, it is conforming, and its boundary matches itself under
the zone's face pairings, which the TetGen mesh does not (GitHub #114, test_18).
Local refinement keeps all three, and respects the resolution limit.
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
    return properties_of(name, _LatticeTri(metric, ops, n))


def properties_of(name, mesh):
    _, ops = primitive_input(name)
    B = primitive_basis(name)
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


def boundary_tetrahedra(name, mesh):
    """Tetrahedra with a face on the zone boundary"""
    B = primitive_basis(name)
    vertices = np.asarray(mesh.vertices) @ B
    tetrahedra = np.asarray(mesh.tetrahedra)
    count = {}
    for t in tetrahedra:
        for f in itertools.combinations(sorted(t), 3):
            count[f] = count.get(f, 0) + 1
    faces = [np.asarray(f) @ B for f in mesh.faces]
    tol = 1e-8 * np.abs(vertices).max()
    planes = []
    for polygon in faces:
        normal = np.linalg.svd(polygon - polygon.mean(axis=0))[2][-1]
        planes.append((normal, normal @ polygon.mean(axis=0)))
    def on_boundary(f):
        return count[f] == 1 and any(all(abs(n @ vertices[k] - d) <= tol for k in f) for n, d in planes)
    return [t for t in tetrahedra if any(on_boundary(f) for f in itertools.combinations(sorted(t), 3))]


@pytest.mark.parametrize("name", CASES)
def test_refined_mesh_properties(name):
    """Refine near a boundary point (boundary triangles split with their partners),
    then repeatedly near the centre with a resolution limit"""
    metric, ops = primitive_input(name)
    B = primitive_basis(name)
    mesh = _LatticeTri(metric, ops, 2)
    vertices = np.asarray(mesh.vertices) @ B
    target = 0.6 * vertices.mean(axis=0) + 0.4 * vertices[np.argmax(np.linalg.norm(vertices, axis=1))]
    for _ in range(3):
        vertices = np.asarray(mesh.vertices) @ B
        near = sorted(boundary_tetrahedra(name, mesh), key=lambda t: np.linalg.norm(vertices[t].mean(axis=0) - target))
        mesh.refine(np.array(near[:6]), 0.0)
    size = np.abs(vertices).max()
    limit = 0.08 * size
    first = len(mesh.vertices)
    for _ in range(4):
        vertices = np.asarray(mesh.vertices) @ B
        centre = vertices.mean(axis=0)
        near = [t for t in np.asarray(mesh.tetrahedra) if np.linalg.norm(vertices[t].mean(axis=0) - centre) < 0.4 * size]
        mesh.refine(np.array(near), limit)
    # no triangle had equivalent longest edges, whose tie only vertex numbering breaks
    assert mesh.self_paired_ties == 0
    result = properties_of(name, mesh)
    assert result["volume"] == pytest.approx(1, abs=1e-9)
    assert result["open_faces"] == 0
    checked, missing = result["boundary"]
    assert checked > 0 and missing == 0, f"{missing} of {checked} boundary images missing"
    # tetrahedra made under the limit have a longest edge above it
    vertices = np.asarray(mesh.vertices) @ B
    new = [t for t in np.asarray(mesh.tetrahedra) if max(t) >= first]
    assert new
    longest = min(max(np.linalg.norm(vertices[a] - vertices[b]) for a, b in itertools.combinations(t, 2)) for t in new)
    assert longest > limit * (1 - 1e-9)


def test_refinement_without_work_changes_nothing():
    metric, ops = primitive_input("hexagonal 6/m")
    mesh = _LatticeTri(metric, ops, 2)
    def state():
        return len(mesh.vertices), {tuple(sorted(t)) for t in np.asarray(mesh.tetrahedra)}
    before = state()
    mesh.refine(np.zeros((0, 4), dtype=int), 0.0)
    assert state() == before
    # nothing is longer than twice this limit
    mesh.refine(np.asarray(mesh.tetrahedra), 1e3)
    assert state() == before
