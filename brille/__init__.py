# Copyright 2020 Greg Tucker
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

"""Python module ``brille``
==============================

This module provides access to the C++ brille library which can be used to
interact with spacegroup and pointgroup symmetry operations, to determine the
first Brillouin zone for a given real space crystallographic lattice, to
find _an_ irreducible polyhedron from the first Brillouin zone and the
pointgroup operations of a spacegroup, to construct polyhedron-filling connected
point networks, and to perform linear-interpolation of user-provided data for
any reciprocal space point at a symmetry-equivalent position within the
connected point network.

.. currentmodule:: brille

.. autosummary::
    :toctree: _generate
"""

# try:
#     from ._brille import *
# except ImportError:
#     print("Error importing compiled brille package")

from ._brille import *
