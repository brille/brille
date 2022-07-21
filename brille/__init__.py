# Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>
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

"""Python module :py:mod:`brille`
=================================

This module provides access to the C++ brille library which can be used to
interact with spacegroup and pointgroup symmetry operations, to determine the
first Brillouin zone for a given real space crystallographic lattice, to
find *an* irreducible polyhedron from the first Brillouin zone and the
pointgroup operations of a spacegroup, to construct polyhedron-filling connected
point networks, and to perform linear-interpolation of user-provided data for
any reciprocal space point at a symmetry-equivalent position within the
connected point network.

.. currentmodule:: brille

.. autosummary::
    :toctree: _generate
"""
from pathlib import Path as _path
from .lattice import wrap_lattice as _wrap_lattice
from .load import load as _load
from . import utils

try:
    from . import plotting
except ModuleNotFoundError:
    # Build servers don't have Matplotlib installed; plotting not tested
    pass


# Find and load the binary brille module
_brille = _load(('_brille',), search=[_path(__file__).parent, _path('.')])

# Store the grid type names for use in, e.g., the plotting routines
__grid_types__ = [f'BZ{x}Q{y}' for x in ('Trellis', 'Mesh', 'Nest') for y in ('cc', 'dc', 'dd')]

# emulate `from _brille import *`
_attrs = [*__grid_types__,
          'AngleUnit',
          'LengthUnit',
          'ApproxConfig',
          'real_space_tolerance',
          'reciprocal_space_tolerance',
          '__version__',
          'version',
          'Basis',
          'Bravais',
          'BrillouinZone',
          'HallSymbol',
          'LPolyhedron',
          'PointSymmetry',
          'Pointgroup',
          'Polyhedron',
          'PrimitiveTransform',
          'RotatesLike',
          'SortingStatus',
          'Spacegroup',
          'Symmetry']
for _attr in _attrs:
    globals()[_attr] = getattr(_brille, _attr)
del _attr, _attrs

# Convert the grid names to their types for use by `isinstance(x, types)`
__grid_types__ = tuple(getattr(_brille, x) for x in __grid_types__)

# wrap _brille.Lattice
Lattice = _wrap_lattice(_brille.Lattice, _brille.Symmetry, _brille.Basis)

# clean-up temporary symbols to avoid accidental exports
del _load, _wrap_lattice, _path
