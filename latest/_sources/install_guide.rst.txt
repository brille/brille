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
If your development environment has `git <https://git-scm.com/>`_, Python ≥ 3.4,
a C++14 compliant compiler, and `CMake <https://cmake.org/>`_ ≥ 3.13,
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
