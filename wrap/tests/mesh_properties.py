"""Properties any mesh of the irreducible Brillouin zone must have.

These checks use only what the compiled module exposes (mesh vertices and
tetrahedra, the zone, the lattice and its symmetry), not how the mesh was built.
They come from the Python prototype of the structured mesh (see the design note)
and serve as acceptance criteria for any mesher:

- the mesh fills the irreducible zone's volume;
- it is conforming: every triangle is shared by two tetrahedra or lies on the
  zone's boundary, so no vertex hangs on another tetrahedron's face or edge;
- its boundary matches itself: a boundary triangle mapped by a symmetry
  operation and a reciprocal lattice translation onto the zone's boundary is a
  boundary triangle too. Otherwise the interpolant is discontinuous across the
  zone faces (GitHub #114).

Coordinates are Cartesian (Å⁻¹), compared to a tolerance relative to the zone.
"""
import itertools
import numpy as np


def mesh_arrays(mesh):
    return np.asarray(mesh.invA, dtype=float), np.asarray(mesh.tetrahedra, dtype=np.int64)


def zone_scale(bz):
    return float(np.abs(np.asarray(bz.ir_vertices_invA)).max())


def ir_planes(bz):
    """Outward face planes (n, d), n·x <= d inside, of the irreducible polyhedron"""
    vertices = np.asarray(bz.ir_vertices_invA, dtype=float)
    centre = vertices.mean(axis=0)
    planes = []
    for face in bz.ir_vertices_per_face:
        polygon = vertices[list(face)]
        normal = np.cross(polygon, np.roll(polygon, -1, axis=0)).sum(axis=0)  # Newell
        normal /= np.linalg.norm(normal)
        d = normal @ polygon.mean(axis=0)
        if normal @ centre > d:
            normal, d = -normal, -d
        planes.append((normal, d))
    return planes


def cartesian_operations(lattice, time_reversal=False):
    """The point group as Cartesian rotations acting on reciprocal vectors"""
    B = np.asarray(lattice.reciprocal_vectors, dtype=float)       # columns
    Bi = np.linalg.inv(B)
    ops = [B @ np.asarray(w, dtype=float).T @ Bi for w in np.asarray(lattice.pointgroup.W)]
    if time_reversal:
        ops += [-r for r in ops]
    unique = []
    for r in ops:
        if not any(np.allclose(r, u, atol=1e-9) for u in unique):
            unique.append(r)
    return unique


def reciprocal_translations(lattice, reach=2):
    """Reciprocal lattice vectors (Cartesian) with coefficients up to `reach`, respecting centring"""
    B = np.asarray(lattice.reciprocal_vectors, dtype=float)
    W = np.asarray(lattice.spacegroup.W)
    w = np.asarray(lattice.spacegroup.w, dtype=float)
    centring = [t for r, t in zip(W, w) if (r == np.eye(3, dtype=int)).all()]
    out = []
    for h in itertools.product(range(-reach, reach + 1), repeat=3):
        h = np.array(h)
        if all(abs(h @ t - round(h @ t)) < 1e-9 for t in centring):
            out.append(B @ h)
    return out


def total_volume(vertices, tetrahedra):
    p = vertices[tetrahedra]
    return float(np.abs(np.einsum('ij,ij->i', p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))).sum() / 6)


def _on_boundary(x, planes, tol):
    return any(abs(n @ x - d) <= tol for n, d in planes)


def _inside(x, planes, tol):
    return all(n @ x - d <= tol for n, d in planes)


def boundary_triangles(vertices, tetrahedra, planes, tol):
    """Triangles used by one tetrahedron that lie in a face plane of the zone"""
    count = {}
    for t in tetrahedra:
        for f in itertools.combinations(sorted(int(i) for i in t), 3):
            count[f] = count.get(f, 0) + 1
    out, other = [], []
    for f, c in count.items():
        if c == 1 and any(all(abs(n @ vertices[i] - d) <= tol for i in f) for n, d in planes):
            out.append(f)
        elif c != 2:
            other.append(f)
    return out, other


def open_faces(vertices, tetrahedra, planes, tol):
    """Triangles used by one tetrahedron but not on the boundary, or by more than two"""
    return len(boundary_triangles(vertices, tetrahedra, planes, tol)[1])


class _PointIndex:
    """Find mesh vertices by position, to a tolerance"""

    def __init__(self, vertices, tol):
        self.tol = tol
        self.cells = {}
        for i, x in enumerate(vertices):
            self.cells.setdefault(tuple(np.floor(x / tol).astype(np.int64)), []).append(i)
        self.vertices = vertices

    def find(self, x):
        base = np.floor(x / self.tol).astype(np.int64)
        for off in itertools.product((-1, 0, 1), repeat=3):
            for i in self.cells.get(tuple(base + off), ()):
                if np.abs(self.vertices[i] - x).max() <= self.tol:
                    return i
        return None


def boundary_mismatches(vertices, tetrahedra, planes, operations, translations, tol):
    """(checked, missing): images of boundary triangles that land on the zone's
    boundary, and how many of those are not boundary triangles of the mesh"""
    tris, _ = boundary_triangles(vertices, tetrahedra, planes, tol)
    present = set(tris)
    index = _PointIndex(vertices, tol)
    checked = missing = 0
    for f in tris:
        p = vertices[list(f)]
        for r in operations:
            q = p @ r.T
            for t in translations:
                img = q + t
                if np.abs(img - p).max() <= tol:
                    continue                                   # the identity map
                if not all(_inside(x, planes, tol) and _on_boundary(x, planes, tol) for x in img):
                    continue
                c = img.mean(axis=0)
                if not any(all(abs(n @ x - d) <= tol for x in img) for n, d in planes) or not _inside(c, planes, tol):
                    continue
                checked += 1
                ids = [index.find(x) for x in img]
                if any(i is None for i in ids) or tuple(sorted(ids)) not in present:
                    missing += 1
    return checked, missing


def check(bz, mesh, lattice, time_reversal=False, relative_tolerance=1e-8):
    """All properties at once, as a dict, for reporting"""
    vertices, tetrahedra = mesh_arrays(mesh)
    tol = relative_tolerance * zone_scale(bz)
    planes = ir_planes(bz)
    checked, missing = boundary_mismatches(vertices, tetrahedra, planes,
                                           cartesian_operations(lattice, time_reversal),
                                           reciprocal_translations(lattice), tol)
    return dict(volume=total_volume(vertices, tetrahedra) / bz.ir_polyhedron.volume,
                open_faces=open_faces(vertices, tetrahedra, planes, tol),
                boundary_checked=checked, boundary_missing=missing)
