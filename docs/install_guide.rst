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

  pip install brille

This will download a precompiled binary of the latest released version of the
module or download its source code and compile the module if a precompiled
binary does not exist for your system.

Note:
  Windows systems must have the Microsoft Visual C++ Redistributable binaries in
  order to load the PyPI distributed :py:mod:`brille` module.
  If you are faced with the error

  .. code-block:: python

    >>> import brille

    ImportError: DLL load failed while importing _brille: The specified module could not be found

  please install the latest x64 Visual C++ redistributable package from
  `Microsoft <https://support.microsoft.com/en-us/help/2977003/the-latest-supported-visual-c-downloads>`_.



source
======
If your development environment has `git <https://git-scm.com/>`_, Python ≥ 3.6,
a C++17 compliant compiler, and `CMake <https://cmake.org/>`_ ≥ 3.13,
then you can build the latest version of brille from source via

.. code-block:: bash

  git clone https://github.com/brille/brille
  cd brille
  python setup.py install

.. role:: bash(code)
  :language: bash
  :class: highlight

Note:
  The CMake build will look for pybind11 and Catch2 header files on your system.
  If they are not present or not a supported version, they will be downloaded
  automatically from their respective Github repositories.

debug
-----
If you encounter the need to debug the :py:mod:`~brille._brille` module, you can
compile and install with debugging symbols via

.. code-block:: bash

   git clone https://github.com/brille/brille
   cd brille
   python setup.py debug_install

which is an alias for ``python setup.py build --debug install``


development
-----------
If you plan to modify the pure Python submodules, e.g.,
:py:mod:`brille.plotting`, you may benefit from
constructing from source using the ``develop`` option.

.. code-block:: bash

    git clone https://github.com/brille/brille
    cd brille
    python setup.py develop

This allows any modifications to the Python source files to be immediately
available to the Python interpreter; where otherwise the ``install`` command
would need to be re-run (which rebuilds the C++ module) each time you would like
to test your modifications.


restricted user access
======================
On some systems the default installation location used by ``pip install`` and
``python setup.py install`` is read-only for standard users.
While one could use an administrator or root account to perform the install in
such a case, a safer alternative is to specify a user-accessible installation
directory via

.. code-block:: bash

  pip install --user brille

or

.. code-block:: bash

  python setup.py --user install


legacy linux systems
====================
If the available compiler and ``pip`` versions are too old and can not be upgraded
you may find that ``pip`` reports that the manylinux2010 versions available
on PyPI are incompatible with your system and building from source may also fail.
This is known to apply to RHEL7 based systems but may affect others as well.

In such a case you can produce a suitable installable package on another system by
replicating the manylinux build system.
For the specific case of RHEL7, starting on a system with `devtoolset-7` installed run

.. code-block:: bash

    scl enable devtoolset-7 bash
    git clone https://github.com/brille/brille.git
    python3 -m pip wheel -w wheelhouse brille

which produces a file like ``brille-0.5.0-cp36-cp36m-linux_x86_64.whl`` that can be copied to the target system.
To install the `brille` package on the target machine one then runs

.. code-block:: bash

    pip install --user brille-0.5.0-cp36-cp36m-linux_x86_64.whl

or similar.
