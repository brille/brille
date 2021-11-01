# Copyright © 2020,2021 Duc Le <duc.le@stfc.ac.uk>
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

"""
Utilities for ``brille``
------------------------

These functions help to construct lattices and grids for brillouin zone
interpolation.

.. code-block:: python

  import numpy as np
  from brille.utils import create_bz, create_grid

  def dispersion(q):
    energies = np.sum(np.cos(np.pi * q), axis=1)
    eigenvectors = np.sin(np.pi * q)
    return energies, eigenvectors
    
  bz = create_bz([4.1, 4.1, 4.1], [90, 90, 90], spacegroup='F m -3 m')
  grid = create_grid(bz, node_volume_fraction=1e-6)
  energies, eigenvectors = dispersion(grid.rlu)
  n_energies = 1     # Only one mode
  n_eigenvectors = 3 # Three values per q
  rotateslike = 0    # See note below
  grid.fill(energies, [n_energies, 0, 0, rotateslike],
            eigenvectors, [n_eigenvectors, 0, 0, rotateslike])

  interp_en, interp_ev = grid.ir_interpolate_at(np.random.rand(1000, 3))

Note
----
The ``rotateslike`` enumeration is given in :py:meth:`~brille._brille.BZMeshQdc.fill`
and describes how the eigenvalues / eigenvectors should be treated on application of a
symmetry operation. The *gamma* option (``rotateslike=3``) should be used for
phonon eigenvectors.


.. currentmodule:: brille.utils

.. autosummary::
    :toctree: _generate
"""

import brille
import numpy as np

def create_bz(*args, is_reciprocal=False, use_primitive=True, search_length=1,
              time_reversal_symmetry=False, wedge_search=True, **kwargs):
    """
    Construct a BrillouinZone object. 

    Parameters
    ----------
    a, b, c : float
        Lattice parameters as separate floating point values
    lens : (3,) :py:class:`numpy.ndarray` or list
        Lattice parameters as a 3-element array or list
    alpha, beta, gamma : float
        Lattice angles in degrees or radians as separate floating point values
        Brille tries to determine if the input is in degrees or radians
        by looking at its magnitude. If the values are all less than PI it
        assumes the angles are in radians otherwise it assumes degrees
    angs : (3,) :py:class:`numpy.ndarray` or list
        Lattice angles in degrees or radians as a 3-element array or list
    lattice_vectors : (3, 3) :py:class:`numpy.ndarray` or list of list
        The lattice vectors as a 3x3 matrix, array or list of list
    spacegroup: str or int
        The spacegroup in either International Tables (Hermann-Mauguin)
        notation or a Hall symbol or an integer Hall number.
    is_reciprocal : bool, keyword-only optional (default: False)
        Whether the lattice parameters or lattice vectors refers to a 
        reciprocal rather than direct lattice. If True, a/b/c/lens should
        be in reciprocal Angstrom, otherwise they should be in Angstrom
    use_primitive : bool, keyword-only optional (default: True)
        Whether the primitive (or conventional) lattice should be used
    search_length : int, keyword-only optional (default: 1)
        An integer to control how-far the vertex-finding algorithm should
        search in τ-index. The default indicates that (1̄1̄1̄), (1̄1̄0), (1̄1̄1),
        (1̄0̄1), ..., (111) are included.
    time_reversal_symmetry : bool, keyword-only optional (default: False)
        Whether to include time-reversal symmetry as an operation to
        determine the irreducible Brillouin zone
    wedge_search : bool, keyword-only optional (default: True)
        If true, return an irreducible first Brillouin zone,
        otherwise just return the first Brillouin zone

    Note
    ----
    Note that the required lattice parameters must be specified as:
        - EITHER ``create_bz(a, b, c, alpha, beta, gamma, spacegroup, ...)``
        - OR     ``create_bz(lens, angs, spacegroup, ...)``
        - OR     ``create_bz(lattice_vectors, spacegroup, ...)``

    E.g. you cannot mix specifing `a`, `b`, `c`, and `angs` etc.
    """
    # Take keyword arguments in preference to positional ones
    a, b, c, alpha, beta, gamma, lens, angs = (kwargs.pop(pname, None)
        for pname in ['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'lens', 'angs'])
    no_lat_kw_s = any([v is None for v in [a, b, c, alpha, beta, gamma]])
    no_lat_kw_v = any([v is None for v in [lens, angs]])
    if no_lat_kw_v and not no_lat_kw_s:
        lens, angs = ([a, b, c], [alpha, beta, gamma])
    lattice_vectors = kwargs.pop('lattice_vectors', None)
    spacegroup = kwargs.pop('spacegroup', None)
    # Parse positional arguments
    spg_id = 0
    if no_lat_kw_s and no_lat_kw_v and lattice_vectors is None:
        if np.shape(args[0]) == ():
            lens, angs = (args[:3], args[3:6])
            spg_id = 6
        elif np.shape(args[0]) == (3,):
            lens, angs = tuple(args[:2])
            spg_id = 2
        elif np.shape(args[0]) == (3,1) or np.shape(args[0]) == (1,3):
            lens, angs = tuple(args[:2])
            lens = np.squeeze(np.array(lens))
            angs = np.squeeze(np.array(angs))
            spg_id = 2
        elif np.shape(args[0]) == (3,3):
            lattice_vectors = args[0]
            spg_id = 1
        else:
            raise ValueError('No lattice parameters or vectors given')
    if spacegroup is None:
        if len(args) > spg_id:
            spacegroup = args[spg_id]
        else:
            raise ValueError('Spacegroup not given')
    if not isinstance(spacegroup, str):
        try:
            spacegroup = int(spacegroup)
        except TypeError as e0:
            e1 = ValueError('Invalid spacegroup input. It must be a string or number')
            e1.__suppress_context__ = True
            e1.__traceback__ = e0.__traceback__
            raise e1

    if is_reciprocal:
        if lattice_vectors is not None:
            lattice = brille.Reciprocal(lattice_vectors, spacegroup)
        else:
            lattice = brille.Reciprocal(lens, angs, spacegroup)
    else:
        if lattice_vectors is not None:
            lattice = brille.Direct(lattice_vectors, spacegroup)
        else:
            lattice = brille.Direct(lens, angs, spacegroup)
        lattice = lattice.star

    try:
        return brille.BrillouinZone(lattice, use_primitive=use_primitive,
                                    search_length=search_length,
                                    time_reversal_symmetry=time_reversal_symmetry,
                                    wedge_search=wedge_search)
    except RuntimeError as e0:
        # We set wedge_search=True by default so add a hint here.
        if 'Failed to find an irreducible Brillouin zone' in str(e0):
            e1 = RuntimeError(str(e0) + ' You can try again with wedge_search=False ' \
                            'to calculate with just the first Brillouin zone')
            e1.__suppress_context__ = True
            e1.__traceback__ = e0.__traceback__
            raise e1
        else:
            raise e0


def create_grid(bz, complex_values=False, complex_vectors=False,
                mesh=False, nest=False, **kwargs):
    """
    Constructs an interpolation grid for a given BrillouinZone object

    Brille provides three different grid implementations:
        - BZTrellisQ: A hybrid Cartesian and tetrahedral grid, with 
          tetrahedral nodes on the BZ surface and cuboids inside. [Default]
        - BZMeshQ: A fully tetrahedral grid with a flat data structure
        - BZNestQ: A fully tetrahedral grid with a nested tree data 
          structure.

    By default a BZTrellisQ grid will be used.

    Parameters
    ----------
    bz : :py:class: `BrillouinZone`
        A BrillouinZone object (required)
    complex_values : bool, optional (default: False)
        Whether the interpolated scalar quantities are complex
    complex_vectors : bool, optional (default: False)
        Whether the interpolated vector quantities are complex
    mesh: bool, optional (default: False)
        Whether to construct a BZMeshQ instead of a BZTrellisQ grid
    nest: bool, optional (default: False)
        Whether to construct a BZNestQ instead of a BZTrellisQ grid

    Note
    ----
    Note that setting both `mesh` and `nest` to True gives an error.


    Additional keyword parameters will be passed to the relevant
    grid constructors.
    
    For ``BZTrellisQ``, these are:

    Parameters
    ----------
    node_volume_fraction : float, optional (default: 1e-5)
        The fractional volume of a tetrahedron in the mesh.
        Smaller numbers will result in better interpolation 
        accuracy at the cost of greater computation time.
    always_triangulate : bool, optional (default: False)
        If set to True, we calculate a bounding polyhedron
        for each point in the grid, and triangulate this into
        tetrahedrons. If False, we set internal points to be 
        cuboid and compute tetrahedrons only for points near
        the surface of the Brillouin Zone.


    For ``BZMeshQ``, these additional parameters are available:

    Parameters
    ----------
    max_size : float, optional (default: -1.0)
        The maximum volume of a tetrahedron in cubic reciprocal
        Angstrom. If set to a negative value, Tetgen will generate
        a tetrahedral mesh without a volume constraint.
    num_levels : int, optional (default: 3)
        The number of layers of triangulation to use.
    max_points : int, optional (default: -1)
        The maximum number of additional mesh points to add to
        improve the mesh quality. Setting this to -1 will allow
        Tetgen to create an unlimited number of additional points.


    For ``BZNestQ``, these additional parameters are available:

    Parameters
    ----------
    max_volume: float
        Maximum volume of a tetrahedron in cubic reciprocal Angstrom.
    number_density: float 
        Number density of points in reciprocal space.
    max_branchings: int, optional (default: 5)
        Maximum number of branchings in the tree structure

    Note
    ----
    Note that one of either the **max_volume** or **number_density**
    parameters must be provided to construct a ``BZNestQ``.
    """
    if not isinstance(bz, brille.BrillouinZone):
        raise ValueError('The `bz` input parameter is not a BrillouinZone object')
    if nest and mesh:
        raise ValueError('Both nest=True and mesh=True is set. Please use one or the other')

    def constructor(grid_type):
        if complex_values and complex_vectors:
            return getattr(brille, grid_type+'cc')
        elif complex_vectors:
            return getattr(brille, grid_type+'dc')
        else:
            return getattr(brille, grid_type+'dd')
 
    if nest:
        if any([v in kwargs for v in ['node_volume_fraction', 'always_triangulate']]):
            raise ValueError('Parameters given are consistent with a trellis grid but nest=True')
        if any([v in kwargs for v in ['max_size', 'num_levels', 'max_points']]):
            raise ValueError('Parameters given are consistent with a mesh grid but mesh=False')
        if 'max_volume' in kwargs:
            return constructor('BZNestQ')(bz, float(kwargs['max_volume']),
                                          kwargs.pop('max_branchings', 5))
        elif 'number_density' in kwargs:
            return constructor('BZNestQ')(bz, int(kwargs['number_density']), 
                                          kwargs.pop('max_branchings', 5))
        else:
            raise ValueError('Neither `max_volume` nor `number_density` provided')
    elif mesh:
        if any([v in kwargs for v in ['node_volume_fraction', 'always_triangulate']]):
            raise ValueError('Parameters given are consistent with a trellis grid but mesh=True')
        if any([v in kwargs for v in ['max_volume', 'max_branchings']]):
            raise ValueError('Parameters given are consistent with a nested grid but nest=False')
        return constructor('BZMeshQ')(bz, float(kwargs.pop('max_size', -1.0)),
                                      int(kwargs.pop('num_levels', 3)),
                                      int(kwargs.pop('max_points', -1)))
    else:
        if any([v in kwargs for v in ['max_volume', 'max_branchings']]):
            raise ValueError('Parameters given are consistent with a nested grid but nest=False')
        if any([v in kwargs for v in ['max_size', 'num_levels', 'max_points']]):
            raise ValueError('Parameters given are consistent with a mesh grid but mesh=False')
        return constructor('BZTrellisQ')(bz, float(kwargs.pop('node_volume_fraction', 1.e-5)),
                                         bool(kwargs.pop('always_triangulate', False)))

