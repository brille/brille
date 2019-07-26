# symbz
A C++ library for symmetry operations with emphasis on the first Brillouin zone.
Wrapped for use in python using [pybind11](https://github.com/pybind/pybind11).

# Dependencies
## TetGen
A slightly modified version of [TetGen](http://tetgen.org) is used to create
refined tetrahedral meshes in the irreducible portion of the first Brillouin
zone.

The modified version is included as part of this repository.

## CGAL
The [Computational Geometry Algorithms Library](https://www.cgal.org/)
(version ≥ 4.13) is used for its quick-location of a point within a tetrahedral
triangulation.

Your linux distribution may provide a packaged version of the CGAL headers.
The current stable Debian (Buster), Fedora (30), and arch distributions provide
compatible header libraries in `libcgal-dev`, `CGAL-devel`, and `cgal` respectively.



Installing CGAL under windows can be tricky so, by
default, the Python build script is set to utilize a CGAL version installed by
[vcpkg](https://github.com/microsoft/vcpkg).
The script looks for an environment variable `VCPKG` which holds the base
installation folder of your vcpkg setup -- you will either need to set a
user or system environment variable with this name and information or set it
on the command line when calling the Python build script.

If you do not have CGAL installed via vcpkg, you can install it via PowerShell

`PS [path to vcpkg]> .\vcpkg.exe install cgal:[arch]`

where `[arch]` is your system architecture triplet.
If the environment variable `VCPKG` is set and you want the 64-bit version of
CGAL, the single-line PowerShell command
`& $Env:VCPKG"\vcpkg.exe" install cgal:x64-windows`
will start the installation for you.

CGAL depends on a number of boost libraries, most of which are installed
automatically as dependencies by vcpkg. At the moment (at least) one boost
dependency, the pointer container library, is missing from the vcpkg
configuration file for CGAL and should be installed via, *e.g.*,
`& $Env:VCPKG"\vcpkg.exe" install boost-ptr-container:x64-windows`.

# Installation
From the root folder of this repository use Python 3 to build and install this
library

`python setup.py install`
