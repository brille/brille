"""Choosing q points to ask brille about.

Only vertices strictly inside the irreducible mesh are safe to test exactly:
vertices on its surface may be moved to a point that is not a vertex (the
mesh does not match itself across symmetry-equivalent faces), or to a point
round-off outside the mesh, which currently crashes brille.
"""
import numpy as np


def interior_vertex_indices(setup):
    """Indices of mesh vertices not on the surface of the tetrahedral mesh."""
    t = np.asarray(setup.grid.tetrahedra)
    faces = np.sort(np.concatenate([t[:, [0, 1, 2]], t[:, [0, 1, 3]], t[:, [0, 2, 3]], t[:, [1, 2, 3]]]), axis=1)
    unique, counts = np.unique(faces, axis=0, return_counts=True)
    surface = np.unique(unique[counts == 1])
    return np.setdiff1d(np.arange(len(setup.vertices)), surface)


def centring_translations(rotations, translations):
    """Non-integer translations of the pure-translation operations (centring vectors)."""
    identity = np.all(np.asarray(rotations).reshape(-1, 9) == np.eye(3, dtype=int).ravel(), axis=1)
    t = np.asarray(translations)[identity] % 1.0
    return t[np.any(np.abs(t - np.round(t)) > 1e-6, axis=1)]


def reciprocal_lattice_shifts(rng, count, max_shift, centring=()):
    """Random reciprocal lattice vectors G, in reciprocal lattice units of the cell.

    For a centred (conventional) cell not every integer triple is a reciprocal
    lattice vector: G must also satisfy ``G·t ∈ ℤ`` for every centring vector t
    (e.g. h + k + l even for an I-centred cell).
    """
    out = np.zeros((0, 3), dtype=int)
    while len(out) < count:
        g = rng.integers(-max_shift, max_shift + 1, size=(2 * count, 3))
        ok = np.ones(len(g), dtype=bool)
        for t in centring:
            gt = g @ t
            ok &= np.abs(gt - np.round(gt)) < 1e-6
        out = np.concatenate([out, g[ok]])
    return out[:count]


def point_group_images(q_rlu, rotations, rng=None, max_shift=0, translations=None):
    """Every point-group image ``Wᵀq`` of each q, optionally shifted by random G.

    ``rotations`` are real-space fractional rotations; reciprocal vectors
    transform with their transposes. Shifts are true reciprocal lattice vectors
    when ``translations`` (of the same operations) reveal a centred cell.
    Returns (Q × |W|, 3).
    """
    w = np.unique(np.asarray(rotations).reshape(-1, 9), axis=0).reshape(-1, 3, 3)
    out = np.einsum("xji,aj->axi", w, np.atleast_2d(q_rlu)).reshape(-1, 3)
    if max_shift:
        centring = () if translations is None else centring_translations(rotations, translations)
        out = out + reciprocal_lattice_shifts(rng, len(out), max_shift, centring)
    return out
