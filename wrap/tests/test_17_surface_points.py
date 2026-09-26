"""Points on the surface of the irreducible zone interpolate.

Moving a point into the irreducible zone can leave it a round-off amount outside
every tetrahedron of the mesh: a symmetry operation or lattice translation maps
it onto a face whose triangulation differs, or the mesh's surface vertices are
themselves rounded. Such points are located in the tetrahedron they are nearest
to, if they are outside it only by round-off (GitHub #121).
"""
import numpy as np
import pytest
from brille import BrillouinZone, BZMeshQdc, Lattice

CASES = {
    # exact values from GitHub #121: two mesh vertices failed to locate
    "triclinic": (((4.023643249400513, 4.90092739265187, 3.2883192254392677),
                   (98.97298894274488, 86.2366290402097, 88.46652897945151)), "-P 1"),
    "fcc": (((5.64,) * 3, (90,) * 3), "-F 4 2 3"),
    "hexagonal": (((3.0, 3.0, 5.0), (90, 90, 120)), "-P 6 2"),
    "monoclinic C": (((5.1, 8.8, 5.2), (90, 104.5, 90)), "-C 2y"),
}


def surface_points(bz, per_face=40, seed=3):
    """Random points on every face of the irreducible polyhedron, then moved by
    random point group operations and lattice translations"""
    rng = np.random.default_rng(seed)
    vertices = np.asarray(bz.ir_vertices)
    points = []
    for face in bz.ir_vertices_per_face:
        polygon = vertices[face]
        for _ in range(per_face):
            # a random point in a random triangle of the fan
            k = rng.integers(1, len(polygon) - 1)
            a, b = sorted(rng.random(2))
            points.append(polygon[0] * a + polygon[k] * (b - a) + polygon[k + 1] * (1 - b))
    points = np.concatenate([np.array(points), vertices])
    rotations = np.asarray(bz.lattice.pointgroup.W).reshape(-1, 3, 3)
    moved = np.einsum("nji,nj->ni", rotations[rng.integers(len(rotations), size=len(points))], points)
    return np.concatenate([points, moved + rng.integers(-2, 3, size=moved.shape)])


@pytest.mark.parametrize("name", CASES)
@pytest.mark.parametrize("threads", (1, 4))
def test_surface_points_interpolate(name, threads):
    lattice, symmetry = CASES[name]
    bz = BrillouinZone(Lattice(lattice, spacegroup=symmetry))
    mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 30)
    mesh.fill(np.ones(len(mesh.rlu)), (1,), mesh.rlu, (0, 3))
    for points in (np.asarray(mesh.rlu), surface_points(bz)):
        values, _ = mesh.ir_interpolate_at(points, threads=threads)
        np.testing.assert_allclose(values, 1, rtol=0, atol=1e-12)
