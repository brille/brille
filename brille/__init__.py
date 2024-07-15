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

from .bound import (
    __version__,
    version,
    AngleUnit,
    LengthUnit,
    Bravais,
    RotatesLike,
    NodeType,
    real_space_tolerance,
    reciprocal_space_tolerance,
    ApproxConfig,
    Basis,
    HallSymbol,
    Spacegroup,
    Pointgroup,
    Symmetry,
    PointSymmetry,
    PrimitiveTransform,
    BrillouinZone,
    Polyhedron,
    LPolyhedron,
    SortingStatus,
    BZTrellisQcc,
    BZTrellisQdc,
    BZTrellisQdd,
    BZMeshQcc,
    BZMeshQdc,
    BZMeshQdd,
    BZNestQcc,
    BZNestQdc,
    BZNestQdd,
    __grid_types__,
)
from .lattice import Lattice
from . import utils
from . import plotting


__all__ = [
    "__version__",
    "version",
    "AngleUnit",
    "LengthUnit",
    "Bravais",
    "RotatesLike",
    "NodeType",
    "real_space_tolerance",
    "reciprocal_space_tolerance",
    "ApproxConfig",
    "Basis",
    "HallSymbol",
    "Spacegroup",
    "Pointgroup",
    "Symmetry",
    "PointSymmetry",
    "PrimitiveTransform",
    "BrillouinZone",
    "Polyhedron",
    "LPolyhedron",
    "SortingStatus",
    "BZTrellisQcc",
    "BZTrellisQdc",
    "BZTrellisQdd",
    "BZMeshQcc",
    "BZMeshQdc",
    "BZMeshQdd",
    "BZNestQcc",
    "BZNestQdc",
    "BZNestQdd",
    "__grid_types__",
    "Lattice",
    "utils",
    "plotting",
]
