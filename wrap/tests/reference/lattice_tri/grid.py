"""The periodic grid of the structured mesh: a lattice's Delaunay triangulation.

Reference implementation for brille's C++ `latticetri::Grid`. Points are exact:
integer coordinates in the grid lattice's basis, scaled by D, so centroids and
(later) midpoints are exact and hashable. A Delaunay triangulation of a lattice is
the Kuhn triangulation of an obtuse superbase (Selling reduction). Where a Selling
parameter is zero the Delaunay cells are larger polytopes; each is split from its
centroid, every face from its own, which is symmetric.
"""
import itertools
import numpy as np

D = 840 * 2 ** 30  # divisible by 1..8, 10, 12, 14 and 30 halvings


def selling(basis):
    """Obtuse superbase (4 vectors summing to zero, pairwise dot <= 0) of the lattice
    whose basis vectors are the rows of `basis` (Cartesian). Returns the superbase as
    Cartesian rows and as integer rows in the original basis."""
    b = [np.array(v, float) for v in basis]
    b.append(-(b[0] + b[1] + b[2]))
    n = [np.array(r) for r in np.eye(3, dtype=int)] + [-np.ones(3, dtype=int)]
    scale = max(np.dot(v, v) for v in b)
    for _ in range(10000):
        pair = next(((i, j) for i in range(4) for j in range(i + 1, 4) if b[i] @ b[j] > 1e-13 * scale), None)
        if pair is None:
            return np.array(b), np.array(n)
        i, j = pair
        k, l = [m for m in range(4) if m not in pair]
        b[k], b[l], b[i] = b[k] + b[i], b[l] + b[i], -b[i]
        n[k], n[l], n[i] = n[k] + n[i], n[l] + n[i], -n[i]
    raise RuntimeError("Selling reduction did not converge")


def kuhn_cell(superbase_int):
    """The 6 Kuhn tetrahedra of one cell, as integer coordinates (lattice basis)."""
    v = superbase_int
    tets = []
    for a, b, c in itertools.permutations(range(3)):
        tets.append((np.zeros(3, int), v[a], v[a] + v[b], v[a] + v[b] + v[c]))
    return tets


class Grid:
    """The lattice with basis rows `basis` (Cartesian), its point group `ops` as
    integer matrices acting on lattice coordinates (columns), and a periodic
    triangulation that is invariant under `ops`."""

    def __init__(self, basis, ops, tie_tolerance=1e-9):
        self.basis = np.array(basis, float)
        self.ops = [np.array(o, int) for o in ops]
        self.metric = self.basis @ self.basis.T
        cart, ints = selling(self.basis)
        self.superbase = cart
        self.selling_parameters = [-(cart[i] @ cart[j]) for i in range(4) for j in range(i + 1, 4)]
        scale = max(np.dot(v, v) for v in cart)
        self.degenerate = min(self.selling_parameters) <= tie_tolerance * scale
        self.cell_tets = kuhn_cell(ints[:3])
        self.tie_tolerance = tie_tolerance

    # exact point keys --------------------------------------------------------
    @staticmethod
    def key(int_coords, denominator=1):
        return tuple(int(x) * (D // denominator) for x in int_coords)

    def cart(self, key):
        return (np.array(key, float) / D) @ self.basis

    def act(self, op, key):
        return tuple(int(x) for x in op @ np.array(key, dtype=object))

    def length2(self, p, q):
        d = np.array([a - b for a, b in zip(p, q)], float) / D
        return float(d @ self.metric @ d)

    # a finite, symmetric patch of the periodic triangulation --------------------
    def patch(self, radius):
        """All cells of the triangulation whose centroid is within `radius`
        (Cartesian) of the origin, as tetrahedra. Whole Delaunay cells are kept or
        dropped, and the ball is invariant under the point group."""
        reach = radius + 2 * max(np.linalg.norm(self.superbase, axis=1))
        # coordinates of points within `reach` are bounded by reach / (smallest singular value)
        span = int(np.ceil(reach / min(np.linalg.svd(self.basis, compute_uv=False)))) + 2
        tets = set()
        for shift in itertools.product(range(-span, span + 1), repeat=3):
            s = np.array(shift)
            for t in self.cell_tets:
                pts = [self.key(p + s) for p in t]
                centroid = sum(self.cart(p) for p in pts) / 4
                if centroid @ centroid <= reach * reach:
                    tets.add(frozenset(pts))
        cells = self.delaunay_cells(tets) if self.degenerate else [[t] for t in tets]
        out = set()
        for members in cells:
            vertices = set().union(*members)
            centroid = sum(self.cart(v) for v in vertices) / len(vertices)
            if centroid @ centroid > radius * radius:
                continue
            out |= set(members) if len(members) == 1 else self._split_cell(members)
        return out

    def delaunay_cells(self, tets):
        """Group Kuhn tetrahedra that share a circumsphere: the Delaunay cells"""
        def circumcentre(t):
            p = [self.cart(k) for k in t]
            a = 2 * np.array([p[1] - p[0], p[2] - p[0], p[3] - p[0]])
            rhs = np.array([p[i] @ p[i] - p[0] @ p[0] for i in (1, 2, 3)])
            return np.linalg.solve(a, rhs)
        scale = max(np.linalg.norm(self.basis, axis=1))
        groups = {}
        for t in tets:
            k = tuple(np.round(circumcentre(t) / (scale * 1e-7)).astype(np.int64))
            groups.setdefault(k, []).append(t)
        return list(groups.values())

    def _split_cell(self, members):
        """Centroid split of the convex cell that is the union of `members`."""
        vertices = sorted(set().union(*members))
        # faces of the cell: triangles of member tetrahedra that are not shared
        count = {}
        for t in members:
            for f in itertools.combinations(sorted(t), 3):
                count[f] = count.get(f, 0) + 1
        boundary = [f for f, c in count.items() if c == 1]
        # merge coplanar boundary triangles into polygons
        pts = {v: self.cart(v) for v in vertices}
        def plane(f):
            a, b, c = (pts[v] for v in f)
            n = np.cross(b - a, c - a)
            n /= np.linalg.norm(n)
            return n, n @ a
        polys = []
        for f in boundary:
            n, d = plane(f)
            for poly in polys:
                if abs(abs(poly['n'] @ n) - 1) < 1e-9 and abs(poly['n'] @ pts[f[0]] - poly['d']) < 1e-9 * (1 + abs(d)):
                    poly['v'] |= set(f)
                    poly['tris'].append(f)
                    break
            else:
                polys.append(dict(n=n, d=d, v=set(f), tris=[f]))
        def centroid(keys):
            keys = list(keys)
            m = len(keys)
            return tuple(sum(k[i] for k in keys) // m for i in range(3)) if all(sum(k[i] for k in keys) % m == 0 for i in range(3)) else None
        c = centroid(vertices)
        assert c is not None, "cell centroid not representable"
        out = set()
        for poly in polys:
            if len(poly['v']) == 3:
                out.add(frozenset(set(poly['v']) | {c}))
                continue
            fc = centroid(poly['v'])
            assert fc is not None, "face centroid not representable"
            # polygon edges: edges of its triangles that occur once
            ec = {}
            for f in poly['tris']:
                for e in itertools.combinations(sorted(f), 2):
                    ec[e] = ec.get(e, 0) + 1
            for e, k in ec.items():
                if k == 1:
                    out.add(frozenset({e[0], e[1], fc, c}))
        return out
