==========
Memory Use
==========

The size of a grid-like object in memory depends on a number of factors, some
of which we can state explicitly and some which can only be estimated in the
general case.

Orthogonal grid
###############

An ortogonal grid in :math:`n` dimensions is comprised of cells which are
defined by their :math:`2^n` vertices.
The :py:mod:`brille` implementation of a three dimenional orthogonal grid is
in :py:class:`BZGridQdd`, :py:class:`BZGridQdc`, and :py:class:`BZGridcc` which
are datatype specialisations of a single C++ class :cpp:class:`BrillouinZoneGrid3`
with double (double), double (complex double),
and complex double (complex double) eigenvalue (eigenvector) element datatypes,
respectively.

The orthogonal grid is guaranteed to envelop the polyhedron used in its
construction and has cell edge lengths determined by either a specified cell
volume or specified number density which is used to calculate the cell volume.
If :math:`V_\text{c}` is the cell volume and :math:`V_\text{p}` the polyhedron
volume, then the cell number density is :math:`N = V_\text{p}/V_\text{c}`.

The number of grid vertices that are within the polyhedron depends on the shape
of the polyhedron, but is similar to the number density if the polyhedron is a
rectangular prism and commensurate with the grid.
In all other cases the number of useful vertices is lower, so :math:`N` is a
useful upper bound.
For each vertex equal-sized data is stored, depending entirely on what the
user provides.
In the case of phonon mode data for a crystal with :math:`n_\text{atom}` atoms
each phonon mode is comprised of one real scalar and :math:`n_\text{atom}`
complex valued three-vectors, and there are :math:`3n_\text{atom}` modes in all.
Thus, for each vertex the grid stores

.. math::

    3 n_\text{atom} + 18 n_\text{atom}^2

8-byte values.

The grid is stored as list of its cell boundaries, each of which is
approximately :math:`N^{1/3}+1` in length, and for each cell at least partly in
the polyhedron the eight corner vertex indices are stored in a structure.
Since all vertices are indexed this meta information totals

.. math::

    N + 3 N^{1/3} + 3

8-byte values.
And therefore, the total memory required for an orthogonal grid storing
phonon data is on the order of

.. math::

    \left[N \left(1 + 3 n_\text{atom} + 18 n_\text{atom}^2\right) + 3 N^{1/3} + 3\right]\times 8\text{bytes}

For a number density of :math:`10^5` a system with ten atoms would therefore
need up to 1.36 GB and one with one hundred atoms would need up to 134 GB.

.. note::

    The actual implementation of :cpp:class:`BZGridQ` handles the meta
    information for in-polyhedron grid vertices differently than described
    here.
    This description is aspirational and would rationalise :cpp:class:`BrillouinZoneGrid3`
    with :cpp:class:`BrillouinZoneTrellis3`.
    Overall memory requirements for the current implementation are similar to
    that described above.

Other grids
###########

For other grids, the per vertex memory usage is similar but the total number of
vertices is not as easy to estimate.
