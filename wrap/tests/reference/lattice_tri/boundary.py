"""The irreducible zone's boundary for the structured mesh, in exact arithmetic.

Reference implementation for brille's C++ `latticetri::Boundary`: the irreducible
polyhedron (first zone ∩ the point group's Dirichlet cone), the maps pairing its
faces, the pairing cells (each face F_j split into the convex pieces F_j ∩ g(F_i)
of positive area), and the special points (cell corners and their images on cell
edges). Coordinates are exact Fractions in the primitive reciprocal lattice basis.
"""
import itertools
from itertools import combinations
from fractions import Fraction as F

import numpy as np


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def solve3(N, d):
    """Exact solution of N x = d (3x3), or None if singular."""
    M = [list(N[i]) + [d[i]] for i in range(3)]
    for c in range(3):
        p = next((r for r in range(c, 3) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        for r in range(3):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return tuple(M[i][3] / M[i][i] for i in range(3))


class Poly:
    """Convex polyhedron as the intersection of halfspaces n.x <= d, exact."""

    def __init__(self, planes):
        self.planes = list(planes)  # (n, d)
        self.verts = {}  # point -> frozenset(plane indices)
        for tri in combinations(range(len(self.planes)), 3):
            x = solve3([self.planes[i][0] for i in tri], [self.planes[i][1] for i in tri])
            if x is None:
                continue
            if all(dot(n, x) <= d for n, d in self.planes):
                self.verts.setdefault(x, set()).update(i for i, (n, d) in enumerate(self.planes) if dot(n, x) == d)
        self.verts = {x: frozenset(s) for x, s in self.verts.items()}

    def adjacent(self, u, v):
        common = self.verts[u] & self.verts[v]
        if len(common) < 2:
            return False
        return not any(common <= s for w, s in self.verts.items() if w != u and w != v)

    def cut(self, n, d):
        k = len(self.planes)
        side = {x: dot(n, x) - d for x in self.verts}
        if all(s <= 0 for s in side.values()):
            return False
        self.planes.append((n, d))
        out = [x for x, s in side.items() if s > 0]
        keep = [x for x, s in side.items() if s < 0]
        new = {}
        for u in out:
            for v in keep:
                if self.adjacent(u, v):
                    t = side[u] / (side[u] - side[v])
                    x = tuple(a + t * (b - a) for a, b in zip(u, v))
                    new.setdefault(x, set()).update(self.verts[u] & self.verts[v])
        verts = {x: s for x, s in self.verts.items() if side[x] <= 0}
        for x, s in verts.items():
            if side[x] == 0:
                verts[x] = s | {k}
        for x, s in new.items():
            verts[x] = frozenset(s | {k})
        self.verts = verts
        return True

    def faces(self):
        return [i for i in range(len(self.planes)) if sum(1 for s in self.verts.values() if i in s) >= 3]


def ir_planes(metric, ops, p=(7, 3, 1)):
    """Half-spaces n·x <= d (Fractions, Λ* coordinates) of the first zone and of
    the point group's Dirichlet cone; and which of them are faces."""
    G = exact_invariant_metric(metric, ops)
    taus = [t for t in itertools.product(range(-2, 3), repeat=3) if any(t)]
    taus.sort(key=lambda t: float(dot(t, [dot(r, t) for r in G])))
    zone = []
    for t in taus:
        n = [dot(r, t) for r in G]
        zone.append((n, dot(t, n) / 2, ('zone', t)))
    M = sum(np.array(o).T @ np.array(o) for o in ops)
    pv = np.array(p)
    cone = []
    for o in ops:
        op = np.array(o) @ pv
        if (op == pv).all():
            continue
        c = M @ (pv - op)
        c = c // np.gcd.reduce(np.abs(c))
        if not any(tuple(c) == tuple(x[0]) for x in cone):
            cone.append((tuple(int(v) for v in c), 'wedge'))
    planes = [(n, d) for n, d, _ in zone] + [([F(-int(v)) for v in c], F(0)) for c, _ in cone]
    tags = [tag for _, _, tag in zone] + [('wedge', c) for c, _ in cone]
    box = [([F(1) if i == j else F(0) for j in range(3)], F(3)) for i in range(3)] + \
          [([F(-1) if i == j else F(0) for j in range(3)], F(3)) for i in range(3)]
    P = Poly(box)
    index = {}                       # position in P.planes -> position in `planes`
    for i, pl in enumerate(planes):
        if P.cut(*pl):
            index[len(P.planes) - 1] = i
    faces = [i for i in P.faces() if i in index]
    assert len(faces) == len(P.faces()), "the box must not bound the zone"
    face_vertices = [[v for v, s in P.verts.items() if i in s] for i in faces]
    return [planes[index[i]] for i in faces], [tags[index[i]] for i in faces], P, face_vertices


def exact_invariant_metric(metric, ops):
    """(1/N) Σ oᵀ G o in exact rationals: invariant under the operations exactly,
    as brille's Lattice now makes its metric. Without it, images of planes miss
    the planes they should coincide with by round-off."""
    G = [[F(float(v)) for v in row] for row in metric]
    out = [[F(0)] * 3 for _ in range(3)]
    for o in ops:
        o = [[int(x) for x in row] for row in o]
        for i in range(3):
            for j in range(3):
                out[i][j] += sum(o[k][i] * G[k][l] * o[l][j] for k in range(3) for l in range(3))
    return [[x / len(ops) for x in row] for row in out]


def side(plane, x):
    n, d = plane
    return dot(n, x) - d


def normalized(plane):
    """Canonical form of the plane {n·x = d} (as a set, orientation ignored)"""
    n, d = plane
    k = next(x for x in n if x != 0)
    return tuple(x / k for x in n), d / k


def image_plane(plane, g, tau):
    """The plane {n·x = d} mapped by x -> g x + tau"""
    n, d = plane
    gi = np.round(np.linalg.inv(np.array(g, float))).astype(int)
    n2 = [sum(n[i] * int(gi[i][k]) for i in range(3)) for k in range(3)]
    return n2, d + dot(n2, tau)


def face_maps(planes, ops, reach):
    """For each face i, the maps (g, tau, j) sending its plane onto face j's plane"""
    face_planes = {normalized(p): i for i, p in enumerate(planes)}
    out = [[] for _ in planes]
    for i, pl in enumerate(planes):
        for g in ops:
            for tau in itertools.product(range(-reach, reach + 1), repeat=3):
                j = face_planes.get(normalized(image_plane(pl, g, tau)))
                if j is not None:
                    out[i].append((g, tau, j))
    return out


def cross(u, v):
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


def clip_convex(poly, window, normal):
    """poly ∩ window, both convex ordered polygons in one plane with this
    (coordinate-space) normal; exact."""
    out = list(poly)
    c = centroid(window)
    for k in range(len(window)):
        a, b = window[k], window[(k + 1) % len(window)]
        n = cross([bi - ai for ai, bi in zip(a, b)], normal)
        if not any(n):
            continue
        d = dot(n, a)
        if dot(n, c) > d:                       # keep the window's side
            n, d = [-x for x in n], -d
        s = [dot(n, v) - d for v in out]
        if all(x <= 0 for x in s):
            continue
        if all(x >= 0 for x in s):
            return []
        out = split_polygon(out, (n, d))[0]     # the side with dot <= d
        out = dedupe(out)
    return out if len(out) >= 3 and exact_rank(out) == 2 else []


def image_points(points, g, tau):
    return [tuple(sum(int(g[i][k]) * x[k] for k in range(3)) + tau[i] for i in range(3)) for x in points]


def pairing_cells(planes, face_vertices, ops, reach=3):
    """For each boundary face F_j, the convex cells F_j ∩ g(F_i) of positive area
    over all maps g (point group op plus translation) and faces F_i, the identity
    on F_j itself excepted. They are the faces of the neighbouring tile(s) across
    F_j, so they partition F_j, and a cell's preimage under g is a cell of F_i:
    cells pair whole, with no closure needed."""
    polys = [order_polygon(fv, np.array([float(c) for c in planes[i][0]])) for i, fv in enumerate(face_vertices)]
    maps = face_maps(planes, ops, reach)
    cells = [[] for _ in planes]
    pairs = [[] for _ in planes]              # (g, tau, i) for each cell of face j
    for i, m in enumerate(maps):
        for g, tau, j in m:
            if i == j and (np.array(g) == np.eye(3, dtype=int)).all() and not any(tau):
                continue
            img = image_points(polys[i], g, tau)
            normal = [F(x) for x in planes[j][0]]
            cell = clip_convex(polys[j], order_polygon(img, np.array([float(c) for c in normal])), normal)
            if cell:
                cells[j].append(order_polygon(cell, np.array([float(c) for c in normal])))
                pairs[j].append((g, tau, i))
    return cells, pairs


def polygon_area2(poly, normal):
    """Twice the area (in coordinate space, times |normal|) of a planar polygon"""
    total = [F(0)] * 3
    for k in range(len(poly)):
        c = cross(poly[k], poly[(k + 1) % len(poly)])
        total = [t + x for t, x in zip(total, c)]
    return abs(dot(total, normal))


def order_polygon(points, normal):
    """Cyclic order of the vertices of a convex polygon in a plane with this normal"""
    P = np.array([[float(c) for c in p] for p in points])
    c = P.mean(axis=0)
    n = normal / np.linalg.norm(normal)
    u = P[0] - c
    u = u - (u @ n) * n
    u /= np.linalg.norm(u)
    w = np.cross(n, u)
    angles = [np.arctan2((q - c) @ w, (q - c) @ u) for q in P]
    return [p for _, p in sorted(zip(angles, points))]


def split_polygon(poly, cut):
    """Split an ordered convex polygon by a plane; exact."""
    s = [side(cut, v) for v in poly]
    if not (any(x > 0 for x in s) and any(x < 0 for x in s)):
        return [poly]
    a, b = [], []
    for i in range(len(poly)):
        p, q = poly[i], poly[(i + 1) % len(poly)]
        sp, sq = s[i], s[(i + 1) % len(poly)]
        if sp <= 0:
            a.append(p)
        if sp >= 0:
            b.append(p)
        if (sp < 0 < sq) or (sq < 0 < sp):
            t = sp / (sp - sq)
            x = tuple(pi + t * (qi - pi) for pi, qi in zip(p, q))
            a.append(x)
            b.append(x)
    return [a, b]


def dedupe(poly):
    """Drop repeated vertices of a cyclic polygon (a line through a corner repeats it)"""
    out = []
    for p in poly:
        if p not in out:
            out.append(p)
    return out


def centroid(points):
    pts = list(points)
    return tuple(sum(p[k] for p in pts) / len(pts) for k in range(3))


def on_segment(y, a, b):
    """y strictly inside segment ab (exact)"""
    d = [bi - ai for ai, bi in zip(a, b)]
    w = [yi - ai for ai, yi in zip(a, y)]
    if any(cross(d, w)):
        return False
    t = dot(w, d) / dot(d, d)
    return 0 < t < 1


def exact_rank(points):
    """Affine rank (0-3) of exact points: 3 unless they are coplanar"""
    p0 = points[0]
    rows = [[a - b for a, b in zip(p, p0)] for p in points[1:]]
    rank = 0
    for col in range(3):
        piv = next((r for r in rows if r[col] != 0), None)
        if piv is None:
            continue
        rows = [[a - (r[col] / piv[col]) * b for a, b in zip(r, piv)] for r in rows if r is not piv]
        rank += 1
    return rank


def special_points(planes, cells, ops):
    """Cell corners, and every image of one (point group operation plus a translation)
    that lies on a cell edge: paired cells then share their vertices along edges."""
    segments = [(cell[k], cell[(k + 1) % len(cell)]) for cs in cells for cell in cs for k in range(len(cell))]
    points = {p for cs in cells for cell in cs for p in cell}
    todo = list(points)
    while todo:
        x = todo.pop()
        for g in ops:
            gx = [sum(int(g[i][k]) * x[k] for k in range(3)) for i in range(3)]
            base = [-int(np.floor(float(c))) for c in gx]
            for tau in itertools.product(*[(b - 1, b, b + 1) for b in base]):
                y = tuple(c + t for c, t in zip(gx, tau))
                if y in points or not all(side(pl, y) <= 0 for pl in planes):
                    continue
                if not any(side(pl, y) == 0 for pl in planes):
                    continue
                if any(on_segment(y, a, b) for a, b in segments):
                    points.add(y)
                    todo.append(y)
    return points


def boundary(metric, ops):
    """(planes, face vertices, cells per face, special points) for a metric and point group"""
    ops = [np.array(o, dtype=int) for o in ops]
    planes, tags, P, face_vertices = ir_planes(np.array(metric, float), ops)
    cells, pairs = pairing_cells(planes, face_vertices, ops)
    return planes, face_vertices, cells, special_points(planes, cells, ops)
