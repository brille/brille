# ![brille]
A C++ library for symmetry operations and linear interpolation within an
irreducible part of the first Brillouin zone.
Wrapped for use in python using [pybind11](https://github.com/pybind/pybind11).

[brille]: https://raw.githubusercontent.com/brille/brille/master/brille.svg

## irreducible polyhedron
When provided with the lattice parameters or basis vectors of a real space spacegroup and its Hall symbol or number as determined by, e.g., [`Spglib`](https://github.com/atztogo/spglib),
`brille` can

- construct the first Brillouin zone
- determine its high symmetry points
- find a irreducible polyhedron and verify that its conforms to the pointgroup symmetry of the spacegroup.

Constructing and irreducible Brillouin zone polyhedron for a face centered cubic lattice can be accomplished with, e.g.,

```python
	import brille

	direct_lattice = brille.Direct((4.96, 4.96, 4.96), (90, 90, 90), 'Fd-3m')
	brillouin_zone = brille.BrillouinZone(direct_lattice.star)
```

## interpolation
Interpolating eigenvalues and eigenvectors across a degenerate point could lead to misidentified equivalent modes.
Since the eigenvectors are distinguishable at all points away from the high-symmetry directions, a hybrid orthogonal/triangulated grid defined in an irreducible part of the Brillouin zone can be used to avoid mode misidentification.

`brille` can construct

- an orthogonal grid guaranteed to contain the first Brillouin zone
- a triangulated set of points filling an irreducible polyhedron or the first Brillouin zone
- or a hybrid orthogonal/triangulated grid where only the surface cells are triangulated to improve point location time while preserving polyhedral conformity.

Of these options the third is most appropriate for the interpolation of models used in the analysis of inelastic neutron scattering data is the hybrid grid.


# Dependencies
## TetGen
A modified version of [TetGen](http://tetgen.org) is used to create
refined tetrahedral meshes in the irreducible portion of the first Brillouin
zone.

The modified version is included as part of this repository.

# Installation
Pre-built wheels for Linux, macOS and Windows are on PyPI:

`python -m pip install brille`

To build from source you need Python 3.11 or later, a C++17 compiler and an
internet connection. From the root folder of this repository run

`python -m pip install .`

The build uses [scikit-build-core](https://scikit-build-core.readthedocs.io)
with [Conan](https://conan.io), which fetches and builds HDF5, HighFive,
pybind11 and Catch2; the first build takes several minutes while HDF5 compiles.
Alternatively, the Python module, C++ library, and [Catch2](https://github.com/catchorg/Catch2)
based tests can be built directly with CMake 3.26 or later; see the
[installation guide](https://brille.github.io/stable/install_guide.html).
