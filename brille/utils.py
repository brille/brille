# Copyright © 2020,2021 Duc Le <duc.le@stfc.ac.uk>
#
# This file is part of brille.
#
# brille is free software: you can redistribute it and/or modify it under the
# terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# brille is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero General Public License
# along with brille. If not, see <https://www.gnu.org/licenses/>.

"""
Utilities for ``brille``
------------------------

These functions help to construct lattices and grids for brillouin zone
interpolation.

.. code-block:: python

  import numpy as np
  from brille.utils import create_bz, create_grid

  def dispersion(q):
    energies = np.sum(np.cos(np.pi * q), axis=1)
    eigenvectors = np.sin(np.pi * q)
    return energies, eigenvectors

  bz = create_bz([4.1, 4.1, 4.1], [90, 90, 90], spacegroup='F m -3 m')
  grid = create_grid(bz, max_size=1e-6)
  energies, eigenvectors = dispersion(grid.rlu)
  n_energies = 1     # Only one mode
  n_eigenvectors = 3 # Three values per q
  rotateslike = 0    # See note below
  grid.fill(energies, [n_energies, 0, 0, rotateslike],
            eigenvectors, [n_eigenvectors, 0, 0, rotateslike])

  interp_en, interp_ev = grid.ir_interpolate_at(np.random.rand(1000, 3))

Note
----
The ``rotateslike`` enumeration is given in :py:meth:`~brille._brille.BZMeshQdc.fill`
and describes how the eigenvalues / eigenvectors should be treated on application of a
symmetry operation. The *gamma* option (``rotateslike=3``) should be used for
phonon eigenvectors.


.. currentmodule:: brille.utils

.. autosummary::
    :toctree: _generate
"""

import numpy as np


def create_bz(
    *args,
    is_reciprocal=False,
    use_primitive=True,
    search_length=1,
    time_reversal_symmetry=False,
    wedge_search=True,
    snap_to_symmetry=True,
    **kwargs,
):
    """
    Construct a BrillouinZone object.

    Parameters
    ----------
    *args
        The lattice and its space group, in one of three forms (see the note
        below): ``a, b, c, alpha, beta, gamma, spacegroup``;
        ``lens, angs, spacegroup``; or ``lattice_vectors, spacegroup``.
    is_reciprocal : bool, keyword-only optional (default: False)
        Whether the lattice parameters or lattice vectors refers to a
        reciprocal rather than direct lattice. If True, a/b/c/lens should
        be in reciprocal Angstrom, otherwise they should be in Angstrom
    use_primitive : bool, keyword-only optional (default: True)
        Whether the primitive (or conventional) lattice should be used
    search_length : int, keyword-only optional (default: 1)
        An integer to control how-far the vertex-finding algorithm should
        search in τ-index. The default indicates that (1̄1̄1̄), (1̄1̄0), (1̄1̄1),
        (1̄0̄1), ..., (111) are included.
    time_reversal_symmetry : bool, keyword-only optional (default: False)
        Whether to include time-reversal symmetry as an operation to
        determine the irreducible Brillouin zone
    wedge_search : bool, keyword-only optional (default: True)
        If true, return an irreducible first Brillouin zone,
        otherwise just return the first Brillouin zone
    snap_to_symmetry: bool, keyword-only optional (default: True)
        Enforces that provided lattice parameters / basis vectors / atom-basis
        positions conform to the provided symmetry operations, if present.

    Other Parameters
    ----------------
    a, b, c : float
        Lattice parameters as separate floating point values
    lens : (3,) :py:class:`numpy.ndarray` or list
        Lattice parameters as a 3-element array or list
    alpha, beta, gamma : float
        Lattice angles in degrees or radians as separate floating point values.
        Brille tries to determine if the input is in degrees or radians
        by looking at its magnitude. If the values are all less than PI it
        assumes the angles are in radians otherwise it assumes degrees
    angs : (3,) :py:class:`numpy.ndarray` or list
        Lattice angles in degrees or radians as a 3-element array or list
    lattice_vectors : (3, 3) :py:class:`numpy.ndarray` or list of list
        The lattice vectors as a 3x3 matrix, array or list of list
    spacegroup : str
        The spacegroup in either International Tables (Hermann-Mauguin)
        notation or a Hall symbol.

    Note
    ----
    The lattice may be given by keyword instead, with the names above; it
    must be given in one of these forms:
        - EITHER ``create_bz(a, b, c, alpha, beta, gamma, spacegroup, ...)``
        - OR     ``create_bz(lens, angs, spacegroup, ...)``
        - OR     ``create_bz(lattice_vectors, spacegroup, ...)``

    E.g. you cannot mix specifing `a`, `b`, `c`, and `angs` etc.
    """
    from brille import Lattice, BrillouinZone

    # Take keyword arguments in preference to positional ones
    a, b, c, alpha, beta, gamma, lens, angs = (
        kwargs.pop(pname, None)
        for pname in ["a", "b", "c", "alpha", "beta", "gamma", "lens", "angs"]
    )
    no_lat_kw_s = any([v is None for v in [a, b, c, alpha, beta, gamma]])
    no_lat_kw_v = any([v is None for v in [lens, angs]])
    if no_lat_kw_v and not no_lat_kw_s:
        lens, angs = ([a, b, c], [alpha, beta, gamma])
    lattice_vectors = kwargs.pop("lattice_vectors", None)
    spacegroup = kwargs.pop("spacegroup", None)
    # Parse positional arguments
    spg_id = 0
    if no_lat_kw_s and no_lat_kw_v and lattice_vectors is None:
        if np.shape(args[0]) == ():
            lens, angs = (args[:3], args[3:6])
            spg_id = 6
        elif np.shape(args[0]) == (3,):
            lens, angs = tuple(args[:2])
            spg_id = 2
        elif np.shape(args[0]) == (3, 1) or np.shape(args[0]) == (1, 3):
            lens, angs = tuple(args[:2])
            lens = np.squeeze(np.array(lens))
            angs = np.squeeze(np.array(angs))
            spg_id = 2
        elif np.shape(args[0]) == (3, 3):
            lattice_vectors = args[0]
            spg_id = 1
        else:
            raise ValueError("No lattice parameters or vectors given")
    if spacegroup is None:
        if len(args) > spg_id:
            spacegroup = args[spg_id]
        else:
            raise ValueError("Spacegroup not given")
    if not isinstance(spacegroup, str):
        e1 = ValueError(
            "Invalid spacegroup input. It must be a string (or valid Spacegroup object)"
        )
        e1.__suppress_context__ = True
        raise e1

    keywords = dict(
        spacegroup=spacegroup,
        real_space=not is_reciprocal,
        snap_to_symmetry=snap_to_symmetry,
    )
    lattice = (
        Lattice((lens, angs), **keywords)
        if lattice_vectors is None
        else Lattice(lattice_vectors, **keywords)
    )

    try:
        return BrillouinZone(
            lattice,
            use_primitive=use_primitive,
            search_length=search_length,
            time_reversal_symmetry=time_reversal_symmetry,
            wedge_search=wedge_search,
        )
    except RuntimeError as e0:
        # We set wedge_search=True by default so add a hint here.
        if "Failed to find an irreducible Brillouin zone" in str(e0):
            e1 = RuntimeError(
                str(e0) + " You can try again with wedge_search=False "
                "to calculate with just the first Brillouin zone"
            )
            e1.__suppress_context__ = True
            e1.__traceback__ = e0.__traceback__
            raise e1
        elif (
            "the dot product between Lattice Vectors requires same or starred lattices"
            in str(e0)
        ):
            e1 = ValueError("Error constructing lattice. Invalid parameters?")
            e1.__suppress_context__ = True
            e1.__traceback__ = e0.__traceback__
            raise e1
        else:
            raise e0


def create_grid(
    bz, complex_values=False, complex_vectors=False, mesh=False, nest=False, trellis=False, **kwargs
):
    """
    Constructs an interpolation grid for a given BrillouinZone object

    Brille provides three different grid implementations:
        - BZMeshQ: A structured tetrahedral mesh, clipped exactly to the zone
          and refinable. [Default]
        - BZTrellisQ: A hybrid Cartesian and tetrahedral grid, with
          tetrahedral nodes on the BZ surface and cuboids inside.
        - BZNestQ: A fully tetrahedral grid with a nested tree data
          structure.

    By default, a BZMeshQ grid is made, unless a ``BZTrellisQ``'s keyword
    arguments are given (as before brille 0.9, when the trellis was the
    default). Without a size, the mesh has ``max_size = 1e-5 / 6``, which gives
    about as many vertices as the trellis's default ``node_volume_fraction =
    1e-5``.

    Parameters
    ----------
    bz : :py:class:`BrillouinZone`
        A BrillouinZone object (required)
    complex_values : bool, optional (default: False)
        Whether the interpolated scalar quantities are complex
    complex_vectors : bool, optional (default: False)
        Whether the interpolated vector quantities are complex
    mesh : bool, optional (default: False)
        Whether to construct a BZMeshQ; it is the default anyway
    nest : bool, optional (default: False)
        Whether to construct a BZNestQ
    trellis : bool, optional (default: False)
        Whether to construct a BZTrellisQ; it is also made when its keyword
        arguments are given

    Other Parameters
    ----------------
    node_volume_fraction : float, optional (default: 1e-5)
        For ``BZTrellisQ``. Despite its name, a volume in cubic reciprocal
        Angstrom, not a fraction: the volume of one cubic node of the trellis,
        which sets its spacing. For a given value, a zone twice the size (in
        each direction) gets eight times the nodes. Smaller numbers give better
        interpolation accuracy at the cost of greater computation time. To size
        the grid by its number of points, use ``bz.ir_polyhedron.volume /
        points``, which gives roughly 1.3 to 2 times ``points`` vertices.
    always_triangulate : bool, optional (default: False)
        For ``BZTrellisQ``. If True, every node is divided into tetrahedra;
        otherwise only the nodes the zone boundary cuts are, and the others stay
        cubes.
    max_size : float, optional (default: -1.0)
        For ``BZMeshQ``. The maximum volume of a grid tetrahedron in cubic
        reciprocal Angstrom, which sets the grid spacing. If not positive, the
        grid is the reciprocal lattice itself, clipped to the zone. Each cell of
        the grid holds six tetrahedra, so a mesh with ``max_size =
        node_volume_fraction / 6`` has about as many vertices as a trellis with
        ``node_volume_fraction`` (up to 1.5 times as many for small grids); and
        ``bz.ir_polyhedron.volume / (6 * points)`` gives roughly 1.5 to 3 times
        ``points`` vertices, the most for small grids.
    num_levels : int, optional (default: 3)
        For ``BZMeshQ``. Unused; kept for compatibility.
    max_points : int, optional (default: -1)
        For ``BZMeshQ``. If positive, the grid is coarsened until its estimated
        number of vertices is at most this.
    max_volume : float
        For ``BZNestQ``. Maximum volume of a tetrahedron in cubic reciprocal
        Angstrom.
    number_density : float
        For ``BZNestQ``. Number density of points in reciprocal space.
    max_branchings : int, optional (default: 5)
        For ``BZNestQ``. Maximum number of branchings in the tree structure.

    Note
    ----
    Setting more than one of `mesh`, `nest` and `trellis` gives an error. Each
    grid's keyword arguments are refused for the other grids, and a ``BZNestQ``
    needs ``nest=True`` and one of **max_volume** or **number_density**.
    """
    from brille import BrillouinZone, _brille

    if not isinstance(bz, BrillouinZone):
        raise ValueError("The `bz` input parameter is not a BrillouinZone object")
    if sum(map(bool, (mesh, nest, trellis))) > 1:
        raise ValueError("Set at most one of mesh=True, nest=True and trellis=True")
    trellis_args = [v for v in ("node_volume_fraction", "always_triangulate") if v in kwargs]
    mesh_args = [v for v in ("max_size", "num_levels", "max_points") if v in kwargs]
    nest_args = [v for v in ("max_volume", "number_density", "max_branchings") if v in kwargs]
    if not (mesh or nest or trellis):
        # the mesh, unless a trellis was asked for by its arguments, as before 0.9
        trellis = bool(trellis_args)
        mesh = not trellis
    chosen = "BZNestQ" if nest else "BZTrellisQ" if trellis else "BZMeshQ"
    for name, given, wanted in (("BZTrellisQ", trellis_args, trellis), ("BZMeshQ", mesh_args, mesh),
                                ("BZNestQ", nest_args, nest)):
        if given and not wanted:
            raise ValueError(f"{', '.join(given)} {'is' if len(given) == 1 else 'are'} for a {name} grid, "
                             f"not the {chosen} being made" + ("; pass nest=True" if name == "BZNestQ" else ""))

    def constructor(grid_type):
        if complex_values and complex_vectors:
            return getattr(_brille, grid_type + "cc")
        elif complex_vectors:
            return getattr(_brille, grid_type + "dc")
        else:
            return getattr(_brille, grid_type + "dd")

    if nest:
        if "max_volume" in kwargs:
            return constructor("BZNestQ")(
                bz, float(kwargs["max_volume"]), kwargs.pop("max_branchings", 5)
            )
        elif "number_density" in kwargs:
            return constructor("BZNestQ")(
                bz, int(kwargs["number_density"]), kwargs.pop("max_branchings", 5)
            )
        else:
            raise ValueError("Neither `max_volume` nor `number_density` provided")
    if trellis:
        return constructor("BZTrellisQ")(
            bz,
            float(kwargs.pop("node_volume_fraction", 1.0e-5)),
            bool(kwargs.pop("always_triangulate", False)),
        )
    # as many vertices, roughly, as the trellis's default: a mesh cell holds six tetrahedra
    default_size = 1.0e-5 / 6 if not mesh_args else -1.0
    return constructor("BZMeshQ")(
        bz,
        float(kwargs.pop("max_size", default_size)),
        int(kwargs.pop("num_levels", 3)),
        int(kwargs.pop("max_points", -1)),
    )


def conventional_to_primitive(lattice, q, values, vectors, tolerance=1e-8):
    """Convert eigen-solutions of a centred conventional cell to its primitive cell

    A grid's eigenvectors describe the atoms of one primitive cell,
    :py:attr:`~brille._brille.Lattice.primitive_basis`, or every atom of the
    conventional cell; convert them to interpolate the primitive cell's modes
    only. A code given a centred
    conventional cell instead returns, at each q, all of its modes: the primitive
    cell's modes at q folded together with those at the other wavevectors the
    larger cell cannot tell from q. This picks out the modes at q, and keeps the
    atoms of the primitive basis.

    Parameters
    ----------
    lattice : :py:class:`~brille._brille.Lattice`
        The conventional lattice, with its full :py:attr:`~brille._brille.Lattice.basis`
        in the order the eigenvectors use.
    q : (Q, 3) array
        The points the eigen-solutions are for, in the conventional cell's
        reciprocal lattice units, as :py:attr:`~brille._brille.BZMeshQdc.rlu`.
    values : (Q, B, ...) array
        The mode values, such as energies, for B = 3 × the conventional cell's atoms
        modes; modes with equal values (in their first element) are degenerate.
    vectors : (Q, B, n, 3) or (Q, B, 3n) array
        The eigenvectors, one 3-vector for each of the conventional cell's n atoms,
        in the cell phase convention: periodic in the conventional cell's reciprocal
        lattice. Their units do not matter.
    tolerance : float, optional
        Relative tolerance for degenerate values.

    Returns
    -------
    values, vectors
        The primitive cell's B / m modes at each q, where m is the number of
        primitive cells in the conventional cell, in the input's order of values,
        with vectors for the atoms of the primitive basis only, in the input's
        layout.

    Note
    ----
    A mode folded from a wavevector p changes by exp(2πi p·d) between an atom and
    its copy d away; the modes at q are those for which p = q, and their primitive
    eigenvector for atom k is the sum over its copies of exp(-2πi q·d) times the
    copy's component, over √m. Degenerate modes can mix the folds, so each
    degenerate group is projected onto the modes at q as a whole.
    """
    import numpy as np

    positions = np.asarray(lattice.basis.positions, dtype=float)
    types = np.asarray(lattice.basis.types)
    primitive = lattice.primitive_basis
    kept_positions = np.asarray(primitive.positions, dtype=float)
    kept_types = np.asarray(primitive.types)
    centring = np.asarray(lattice.centring_vectors, dtype=float)
    m, n_primitive, n_conventional = len(centring), len(kept_positions), len(positions)
    if n_conventional != m * n_primitive:
        raise ValueError(f"the basis has {n_conventional} atoms, not {m} copies of the {n_primitive}-atom primitive cell")
    # each atom's primitive atom, and the actual vector to it from that atom
    owner = np.empty(n_conventional, dtype=int)
    offset = np.empty((n_conventional, 3))
    for j, (p, t) in enumerate(zip(positions, types)):
        for k in np.flatnonzero(kept_types == t):
            d = p - kept_positions[k]
            r = d[None, :] - centring
            if np.any(np.all(np.abs(r - np.round(r)) < 1e-8, axis=1)):
                owner[j], offset[j] = k, d
                break
        else:
            raise ValueError(f"atom {j} is not a centring copy of an atom of the primitive basis")

    q = np.atleast_2d(np.asarray(q, dtype=float))
    values = np.asarray(values)
    vectors = np.asarray(vectors)
    shape = vectors.shape
    vectors = vectors.reshape(shape[0], shape[1], n_conventional, 3)
    n_q, n_modes = vectors.shape[:2]
    if n_modes % m:
        raise ValueError(f"{n_modes} modes are not {m} folds of the primitive cell's")
    phase = np.exp(-2j * np.pi * q @ offset.T)                    # (Q, n)
    fold = np.zeros((n_primitive, n_conventional))
    fold[owner, np.arange(n_conventional)] = 1.0
    out_values = np.empty((n_q, n_modes // m) + values.shape[2:], dtype=values.dtype)
    out_vectors = np.empty((n_q, n_modes // m, n_primitive, 3), dtype=complex)
    for i in range(n_q):
        key = values[i].reshape(n_modes, -1)[:, 0].real
        order = np.argsort(key, kind="stable")
        kept = 0
        start = 0
        while start < n_modes:
            stop = start + 1
            while stop < n_modes and abs(key[order[stop]] - key[order[start]]) <= tolerance * max(1.0, abs(key[order[start]])):
                stop += 1
            group = order[start:stop]
            modes = vectors[i, group]                                   # (g, n, 3)
            # project onto the modes at q: (P e)_j = exp(2πi q.d_j) mean over the copies of
            # the atom of exp(-2πi q.d) e. P maps a degenerate group into itself, and in the
            # group's own coordinates it is a Hermitian projector, whatever the units.
            summed = np.einsum("kj,j,gjd->gkd", fold, phase[i], modes) / m             # (g, N, 3)
            projected = np.einsum("j,gjd->gjd", np.conj(phase[i]), summed[:, owner])   # (g, n, 3)
            flat = modes.reshape(len(group), -1).T
            c = np.linalg.lstsq(flat, projected.reshape(len(group), -1).T, rcond=None)[0]
            w, c = np.linalg.eigh((c + c.conj().T) / 2)
            at_q = c[:, w > 0.5]
            if kept + at_q.shape[1] > n_modes // m:
                break
            for r in range(at_q.shape[1]):
                mode = np.einsum("g,gjd->jd", at_q[:, r], modes)      # a mode at q, in the input's units
                out_values[i, kept] = values[i, group[0]]
                out_vectors[i, kept] = np.einsum("kj,j,jd->kd", fold, phase[i], mode) / np.sqrt(m)
                kept += 1
            start = stop
        if kept != n_modes // m:
            raise ValueError(f"found {kept} of the primitive cell's {n_modes // m} modes at q = {q[i]}: are the "
                             "eigenvectors in the cell phase convention, for the atoms of lattice.basis in order?")
    if len(shape) == 3:
        out_vectors = out_vectors.reshape(n_q, n_modes // m, 3 * n_primitive)
    return out_values, out_vectors
