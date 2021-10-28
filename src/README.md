A C++ library for symmetry operations and linear interpolation within an
irreducible part of the first Brillouin zone.

## irreducible polyhedron
When provided with the lattice parameters or basis vectors of a real space lattice and its spacegroup symmetry operators (explicitly, via a Hall symbol, or as CIF xyz strings) `brille` can

- construct the first Brillouin zone
- determine its high symmetry points
- find a irreducible polyhedron and verify that its conforms to the pointgroup
  symmetry of the spacegroup.

## interpolation
Interpolating eigenvalues and eigenvectors across a degenerate point could lead to misidentified equivalent modes. Since the eigenvectors are distinguishable at all points away from the high-symmetry directions, a hybrid orthogonal/triangulated grid defined in an irreducible part of the Brillouin zone can be used to avoid mode misidentification.

`brille` can construct

- a triangulated set of points filling an irreducible polyhedron or the first Brillouin zone
- or a hybrid orthogonal/triangulated grid where only the surface cells are triangulated to improve point location time while preserving polyhedral conformity.

Of these options that most appropriate for the interpolation of models used in the analysis of inelastic neutron scattering data is the hybrid grid.


# Dependencies
## TetGen
A modified version of [TetGen](http://tetgen.org) is used to create refined tetrahedral meshes in the irreducible portion of the first Brillouin zone. The modified version is included as part of the `brille` library header.

# Installation
Use `CMake` to configure and build the `brille` single header `brille.h` and its shared or static associated library file.
From the repository root directory, run
```bash
  cmake -S . -B build
  cmake --build build --target library
```

The library components can then be installed into the default location by executing
```bash
  cmake --install build
```
as a user with elevated privileges. Alternatively, the installation directory can be changed to a location writable by your user, e.g.,
```bash
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX="~/.local"
  cmake --install build
```

# Documentation
The C++ library documentation can be built if your system has a working Doxygen installation. By default the documentation target is turned off so you must turn it on by
```bash
  cmake -S . -B build -DUSE_DOXYGEN=ON
```
after which it the `CMake` target 'docs' can be used to build the HTML documentation.

```bash
  cmake --build build --target docs
```
