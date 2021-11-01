==========
Memory Use
==========

The size of a grid-like object in memory depends on a number of factors, some
of which we can state explicitly and some which can only be estimated in the
general case.

Orthogonal grid
###############

An orthogonal grid in :math:`n` dimensions is comprised of cells which are
defined by their :math:`2^n` vertices.
While :py:mod:`brille` does not implement a basic three dimensional orthogonal
grid, it represents a useful starting point for memory use calculations.

An orthogonal grid would need to envelop the polyhedron used in its
construction and would have cell edge lengths determined by either a specified
cell volume or specified number density.
If :math:`V_\text{c}` is the cell volume and :math:`V_\text{p}` the polyhedron
volume, then the number of cells is :math:`N = V_\text{p}/V_\text{c}`.

The number of grid vertices within the polyhedron would depend on the shape
of the polyhedron, but should be similar to the number density times the
polyhedron volume in all cases.
As with the :py:mod:`brille` grids, for each vertex equal-sized data is stored,
depending entirely on what the user provides.
In the case of phonon mode data for a crystal with :math:`n_\text{atom}` atoms
each phonon mode is comprised of one real scalar and :math:`n_\text{atom}`
complex valued three-vectors, and there are :math:`3n_\text{atom}` modes in all.
Thus, for each vertex the grid stores
:math:`3 n_\text{atom} + 18 n_\text{atom}^2` 8-byte values.

The grid itself would be stored as list of its cell boundaries, each of which is
approximately :math:`N^{1/3}+1` in length, and for each cell at least partly in
the polyhedron the eight corner vertex indices would be stored stored in a
structure as well.
Since all vertices are indexed this meta information totals
:math:`N + 3 N^{1/3} + 3` 8-byte values.
And therefore, the total memory required for an orthogonal grid storing
phonon data is on the order of

.. math::
    \left[N \left(1 + 3 n_\text{atom} + 18 n_\text{atom}^2\right) + 3 N^{1/3} + 3\right]\times 8\text{bytes}
    :label: memuse

For a grid with :math:`10^5` cells, a system with ten atoms would therefore
need up to 1.36 GB and one with one hundred atoms would need up to 134 GB.

Other grids
###########

For other grids, the per vertex memory usage is similar but the total number of
vertices is not as easy to estimate. The expression in :eq:`memuse` should be
an appropriate order-of-magnitude estimate of memory use.
