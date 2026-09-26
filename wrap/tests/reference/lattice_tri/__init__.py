"""Reference implementation of the structured mesh's building blocks (see the
design note), in exact arithmetic. It is slow, and only used to check the C++.

- grid.py: Selling reduction, Kuhn tetrahedra, and the symmetric split of
  degenerate Delaunay cells, i.e. the periodic grid.
"""
