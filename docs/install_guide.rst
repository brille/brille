.. _install_guide:

==================
Installation Guide
==================
This module relies heavily on C++ wrapped using `pybind11 <https://github.com/pybind/pybind11>`_.
You can obtain a working copy through the `Python Package Index <https://pypi.org/>`_ using ``pip``
or by compiling the source after downloading it from the `brille Github repository <https://github.com/brille/brille>`_.

pip
===
With an up-to-date version of ``pip`` you can install brille via

.. code-block:: bash

  python -m pip install brille

This downloads a precompiled wheel of the latest release, or downloads the
source code and compiles it if no wheel matches your system (see
`building from source`_). Wheels are published for

* Linux on x86-64 with glibc 2.28 or later (manylinux_2_28), or with musl libc
  (musllinux_1_2),
* macOS 14 or later, on Apple silicon and Intel,
* Windows on x86-64.

The plotting functions in :py:mod:`brille.plotting` need matplotlib, which you
can install alongside via ``python -m pip install "brille[plotting]"``.

Note:
  Windows systems must have the Microsoft Visual C++ Redistributable binaries in
  order to load the PyPI distributed :py:mod:`brille` module.
  If you are faced with the error

  .. code-block:: python

    >>> import brille

    ImportError: DLL load failed while importing _brille: The specified module could not be found

  please install the latest x64 Visual C++ redistributable package from
  `Microsoft <https://support.microsoft.com/en-us/help/2977003/the-latest-supported-visual-c-downloads>`_.


building from source
====================
You need Python 3.11 or later, a C++17 compiler, and an internet connection.
Then build and install the latest version of brille with

.. code-block:: bash

  git clone https://github.com/brille/brille
  cd brille
  python -m pip install .

The build uses `scikit-build-core <https://scikit-build-core.readthedocs.io>`_
with `Conan <https://conan.io>`_, which fetches and builds HDF5, HighFive,
pybind11 and Catch2 from ConanCenter; ``pip`` installs both into an isolated
build environment, so keep build isolation on (don't pass
``--no-build-isolation``). The first build takes several minutes while HDF5
compiles; later builds reuse Conan's cache.

debugging symbols
-----------------
The build is always an optimised release build, and the module is stripped.
On Linux and macOS you can keep debugging symbols, for use with a debugger or a
native profiler such as ``py-spy --native``, with

.. code-block:: bash

   python -m pip install . -C cmake.define.CMAKE_CXX_FLAGS=-g -C cmake.define.CMAKE_STRIP=/bin/true


development
-----------
If you plan to modify the pure Python submodules, e.g.,
:py:mod:`brille.plotting`, install in editable mode:

.. code-block:: bash

    python -m pip install -e .

Changes to the Python source files are then available immediately. Changes to
the C++ source still need the command to be run again, which rebuilds the
module.

CMake and the C++ tests
-----------------------
The Python module, the C++ library and the `Catch2 <https://github.com/catchorg/Catch2>`_
based tests can also be built directly with `CMake <https://cmake.org/>`_ 3.26 or later.
CMake runs Conan itself, so install it into the Python environment first:

.. code-block:: bash

    python -m pip install conan setuptools_scm numpy
    cmake -S . -B build -D CMAKE_BUILD_TYPE=Release -D Python3_EXECUTABLE=$(which python)
    cmake --build build --config Release -j
    ctest --test-dir build -C Release


restricted user access
======================
On some systems the default installation location used by ``pip install`` is
read-only for standard users.
While one could use an administrator or root account to perform the install in
such a case, a safer alternative is a virtual environment, or a user-accessible
installation directory via

.. code-block:: bash

  python -m pip install --user brille
