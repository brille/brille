"""Refining a BZMeshQdc in place: refinement_points shows what refine will add,
refine appends vertices (and their data) without disturbing existing ones, and
the mesh keeps the properties of test_18, stops at a resolution limit, and
survives an HDF5 round trip.
"""
import itertools

import numpy as np
import pytest
from brille import BrillouinZone, BZMeshQdc, Lattice
from mesh_properties import check

LATTICES = {
    "hexagonal 6/m": (((3.0, 3.0, 5.0), (90, 90, 120)), "-P 6"),
    "triclinic -P 1": (((4.02, 4.90, 3.29), (98.97, 86.24, 88.47)), "-P 1"),
    "fcc": (((5.64,) * 3, (90,) * 3), "-F 4 2 3"),
}


def zone(name):
    lengths_angles, symmetry = LATTICES[name]
    lattice = Lattice(lengths_angles, spacegroup=symmetry)
    return lattice, BrillouinZone(lattice)


def smooth(q):
    """A periodic function to interpolate"""
    q = np.asarray(q)
    return np.cos(2 * np.pi * q[:, 0]) + np.cos(2 * np.pi * q[:, 1]) + np.cos(2 * np.pi * q[:, 2])


def filled(bz, divisor=200):
    mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / divisor)
    rlu = np.asarray(mesh.rlu)
    mesh.fill(smooth(rlu)[:, None], (1,), np.zeros((len(rlu), 1), dtype=complex), (1,))
    return mesh


def data_for(points):
    return smooth(points)[:, None], np.zeros((len(points), 1), dtype=complex)


def longest_edges(mesh, tetrahedra):
    vertices = np.asarray(mesh.invA)
    return np.array([max(np.linalg.norm(vertices[a] - vertices[b]) for a, b in itertools.combinations(t, 2)) for t in tetrahedra])


def test_refinement_points_are_what_refine_adds():
    _, bz = zone("hexagonal 6/m")
    mesh = filled(bz)
    before = np.asarray(mesh.rlu).copy()
    tetrahedra = np.asarray(mesh.tetrahedra).copy()
    where = np.arange(0, len(tetrahedra), 7)
    planned = mesh.refinement_points(where)
    assert len(planned) > 0
    # the dry run changes nothing
    assert np.array_equal(np.asarray(mesh.rlu), before)
    assert np.array_equal(np.asarray(mesh.tetrahedra), tetrahedra)
    added = mesh.refine(where, *data_for(planned))
    assert np.array_equal(added, planned)
    # new vertices come after the existing ones, which keep their indices
    rlu = np.asarray(mesh.rlu)
    assert np.array_equal(rlu[:len(before)], before)
    assert np.array_equal(rlu[len(before):], planned)


def test_refine_keeps_and_extends_the_data():
    _, bz = zone("hexagonal 6/m")
    mesh = filled(bz)
    where = np.arange(0, len(mesh.tetrahedra), 5)
    mesh.refine(where, *data_for(mesh.refinement_points(where)))
    rlu = np.asarray(mesh.rlu)
    values = np.asarray(mesh.ir_interpolate_at(rlu, threads=1)[0]).ravel()
    assert values == pytest.approx(smooth(rlu), abs=1e-12)


def test_refine_refuses_data_that_does_not_fit():
    _, bz = zone("hexagonal 6/m")
    mesh = filled(bz)
    where = np.arange(0, len(mesh.tetrahedra), 7)
    planned = mesh.refinement_points(where)
    count = len(mesh.rlu)
    with pytest.raises(RuntimeError, match="adds"):
        mesh.refine(where, *data_for(planned[:-1]))
    with pytest.raises(RuntimeError, match="holds data"):
        mesh.refine(where)
    assert len(mesh.rlu) == count
    empty = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 200)
    with pytest.raises(RuntimeError, match="no data"):
        empty.refine(where, *data_for(planned))


@pytest.mark.parametrize("where", [np.array([10 ** 9]), np.zeros(3, dtype=bool), np.array([[0, 1]])])
def test_where_must_name_tetrahedra(where):
    _, bz = zone("fcc")
    mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 100)
    with pytest.raises((IndexError, ValueError)):
        mesh.refinement_points(where)


@pytest.mark.parametrize("name", LATTICES)
def test_refined_mesh_properties(name):
    """Partial refinement, repeated, keeps the mesh filling the zone, conforming,
    and matching itself across paired zone faces"""
    lattice, bz = zone(name)
    mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 100)
    for step in range(3):
        mask = np.zeros(len(mesh.tetrahedra), dtype=bool)
        mask[step::4] = True
        mesh.refine(mask)
    result = check(bz, mesh, lattice)
    assert result["volume"] == pytest.approx(1, abs=1e-9)
    assert result["open_faces"] == 0
    assert result["boundary_missing"] == 0, f"{result['boundary_missing']} of {result['boundary_checked']} boundary images missing"


def test_refinement_stops_at_the_resolution_limit():
    _, bz = zone("hexagonal 6/m")
    mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 200)
    first = len(mesh.rlu)
    limit = 0.1   # resolution / points_per_resolution, inverse Angstrom
    for _ in range(20):
        if not len(mesh.refine(resolution=2 * limit, points_per_resolution=2)):
            break
    else:
        pytest.fail("refinement did not stop")
    tetrahedra = np.asarray(mesh.tetrahedra)
    made = tetrahedra[tetrahedra.max(axis=1) >= first]
    longest = longest_edges(mesh, made)
    # nothing is split below the limit, and nothing left above twice it
    assert longest.min() > limit
    assert longest.max() <= 2 * limit * (1 + 1e-12)


def test_uniform_refinement_improves_interpolation():
    _, bz = zone("fcc")
    mesh = filled(bz, 100)
    q = np.random.default_rng(1).uniform(-1, 1, (5000, 3))
    exact = smooth(q)

    def error():
        return np.sqrt(np.mean((np.asarray(mesh.ir_interpolate_at(q, threads=1)[0]).ravel() - exact) ** 2))

    before = error()
    mesh.refine(None, *data_for(mesh.refinement_points()))
    assert error() < 0.6 * before


@pytest.mark.parametrize("name", LATTICES)
def test_repeated_uniform_refinement(name):
    """Tetrahedra and boundary triangles break longest-edge ties alike; when they
    didn't, the hexagonal mesh's third uniform refinement propagated forever"""
    lattice, bz = zone(name)
    mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 200)
    for _ in range(3):
        assert len(mesh.refine(None))
    result = check(bz, mesh, lattice)
    assert result["volume"] == pytest.approx(1, abs=1e-9)
    assert result["open_faces"] == 0
    assert result["boundary_missing"] == 0, f"{result['boundary_missing']} of {result['boundary_checked']} boundary images missing"


def test_refinement_survives_saving(tmp_path):
    _, bz = zone("triclinic -P 1")
    mesh = filled(bz)
    where = np.arange(0, len(mesh.tetrahedra), 9)
    mesh.refine(where, *data_for(mesh.refinement_points(where)))
    path = tmp_path / "mesh.h5"
    mesh.to_file(str(path))
    loaded = BZMeshQdc.from_file(str(path))
    assert loaded.refinable
    assert np.array_equal(np.asarray(loaded.rlu), np.asarray(mesh.rlu))
    assert np.array_equal(np.asarray(loaded.tetrahedra), np.asarray(mesh.tetrahedra))
    again = np.arange(0, len(mesh.tetrahedra), 13)
    planned = mesh.refinement_points(again)
    assert np.array_equal(loaded.refinement_points(again), planned)
    assert np.array_equal(loaded.refine(again, *data_for(planned)), mesh.refine(again, *data_for(planned)))


def test_released_triangulation_is_rebuilt_alike():
    _, bz = zone("hexagonal 6/m")
    kept = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 200)
    released = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 200)
    assert not kept.holds_triangulation
    for step in (7, 5):
        where = np.arange(0, len(kept.tetrahedra), step)
        assert np.array_equal(kept.refine(where), released.refine(where))
        released.release_triangulation()
        assert kept.holds_triangulation and not released.holds_triangulation
    where = np.arange(0, len(kept.tetrahedra), 3)
    assert np.array_equal(kept.refinement_points(where), released.refinement_points(where))
    assert np.array_equal(np.asarray(kept.tetrahedra), np.asarray(released.tetrahedra))
