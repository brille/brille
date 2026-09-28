"""
    pybind11 module :py:mod:`brille._brille`
    ----------------------------------------
    This module provides the interface to the C++ library.

    All of the symbols defined within :py:mod:`brille._brille` are imported by
    :py:mod:`brille` to make using them easier.
    If in doubt, the interfaced classes can be accessed via their submodule
    syntax.

    .. code-block:: python

      from brille._brille import Direct, BrillouinZone
      from brille.plotting import plot as bplot

      direct_lattice = Direct((3.95, 3.95, 3.95, 12.9), (90, 90, 90), 'I4/mmm')
      brillouin_zone = BrillouinZone(direct_lattice.star)

      bplot(brillouin_zone)

    .. currentmodule:: brille._brille

    .. autosummary::
      :toctree: _generate

  
"""
from __future__ import annotations
import numpy
import pybind11_stubgen.typing_ext
import typing
__all__: list[str] = ['AngleUnit', 'ApproxConfig', 'BZMeshQcc', 'BZMeshQdc', 'BZMeshQdd', 'BZNestQcc', 'BZNestQdc', 'BZNestQdd', 'BZTrellisQcc', 'BZTrellisQdc', 'BZTrellisQdd', 'Basis', 'Bravais', 'BrillouinZone', 'HallSymbol', 'LPolyhedron', 'Lattice', 'LengthUnit', 'NearSymmetryWarning', 'NodeType', 'PointSymmetry', 'Pointgroup', 'Polyhedron', 'PrimitiveTransform', 'RotatesLike', 'SortingStatus', 'Spacegroup', 'Symmetry', 'build_datetime', 'build_hostname', 'emit', 'emit_datetime', 'git_branch', 'git_revision', 'real_space_tolerance', 'reciprocal_space_tolerance', 'version']
class NearSymmetryWarning(UserWarning):
    """
    The lattice is close to one with more symmetry, so its Brillouin zone has features much smaller than itself.
    """
class AngleUnit:
    """
      The units of a number representing an angle.
    
      >>> from brille import AngleUnit
      >>> r = AngleUnit.radian
      
    
    Members:
    
      not_provided : unknown, will be inferred
    
      radian : radian
    
      degree : degree
    
      pi : radian divided by pi
    """
    __members__: typing.ClassVar[dict[str, AngleUnit]]
    degree: typing.ClassVar[AngleUnit]
    not_provided: typing.ClassVar[AngleUnit]
    pi: typing.ClassVar[AngleUnit]
    radian: typing.ClassVar[AngleUnit]
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class LengthUnit:
    """
      The units of a number representing a length.
    
      >>> from brille import LengthUnit
      >>> a = LengthUnit.angstrom
      
    
    Members:
    
      none : A default value, the used value may be inferred
    
      angstrom : 10\\ :sup:`-10` meter
    
      inverse_angstrom : 1 / angstrom
    
      real_lattice : fractional coordinates of real lattice
    
      reciprocal_lattice : fractional coordinates of reciprocal lattice
    """
    __members__: typing.ClassVar[dict[str, LengthUnit]]
    angstrom: typing.ClassVar[LengthUnit]
    inverse_angstrom: typing.ClassVar[LengthUnit]
    none: typing.ClassVar[LengthUnit]
    real_lattice: typing.ClassVar[LengthUnit]
    reciprocal_lattice: typing.ClassVar[LengthUnit]
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class NodeType:
    """
      An enumeration to differentiate between Node types of the :py:class:`brille._brille.BZTrellisQdd`,
      :py:class:`brille._brille.BZTrellisQdc`, and :py:class:`brille._brille.BZTrellisQcc`, classes.
    
      The return type of the methods, e.g.,
      :py:meth:`~brille._brille.BZTrellisQdc.node_at_type`,
      :py:meth:`~brille._brille.BZTrellisQdc.node_containing_type`, and
      :py:meth:`~brille._brille.BZTrellisQdc.all_node_types`.
    
      >>> from brill import NodeType
      >>> nt = NodeType.polygon
      
    
    Members:
    
      assumed_null : Indicates a node which should be null but no explicit check was performed
    
      found_null : A node which an explicit check found to be null
    
      null : A null node
    
      cube : A cube shaped node
    
      polygon : A convex polygon shaped node
    """
    __members__: typing.ClassVar[dict[str, NodeType]]
    assumed_null: typing.ClassVar[NodeType]
    cube: typing.ClassVar[NodeType]
    found_null: typing.ClassVar[NodeType]
    null: typing.ClassVar[NodeType]
    polygon: typing.ClassVar[NodeType]
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class ApproxConfig:
    """
    A local set of approximate floating point comparison values
    """
    def __init__(self, digits: int = 1000, real_space_tolerance: float = 1e-10, reciprocal_space_tolerance: float = 1e-10) -> None:
        ...
    @property
    def digits(self) -> int:
        """
        A multiplier on the machine epsilon, below which two floating point numbers are approximately the same.
        """
    @digits.setter
    def digits(self, arg1: int) -> int:
        ...
    @property
    def real_space_tolerance(self) -> float:
        """
        An absolute real space floating point tolerance in angstrom
        """
    @real_space_tolerance.setter
    def real_space_tolerance(self, arg1: float) -> float:
        ...
    @property
    def reciprocal_space_tolerance(self) -> float:
        """
        An absolute reciprocal space floating point tolerance in inverse angstrom
        """
    @reciprocal_space_tolerance.setter
    def reciprocal_space_tolerance(self, arg1: float) -> float:
        ...
class Bravais:
    """
      A Bravais letter indicating the centering of a lattice
    
      When the unit cell does not reflect the symmetry of the lattice, it is usual
      to refer to a 'conventional' crystallographic basis,
      :math:`(\\mathbf{a}_s\\,\\mathbf{b}_s\\,\\mathbf{c}_s)`, instead of
      a primitive basis, :math:`(\\mathbf{a}_p\\,\\mathbf{b}_p\\,\\mathbf{c}_p)`.
      Such a conventional basis has 'extra' lattice points added at the centre of
      the unit cell, the centre of a face, or the centre of three faces.
      The 'extra' nodes in the conventional basis are displaced from the origin of
      the unit cell by 'centring vectors'. As with any space-spanning basis, any
      whole-number linear combination of the conventional basis vectors is a lattice
      point but in addition there exist linear combinations
      :math:`x\\mathbf{a}_s+y\\mathbf{b}_s+z\\mathbf{c}_s` with at least
      two fractional coefficients :math:`(x,y,z)` that are lattice points as well.
    
      Each conventional basis is ascribed a Bravais letter, which forms part of the
      Hermann-Mauguin symbol of a space group.
      A subset of the 10 possible Bravais letters is used herein:
    
      | Bravais letter | Centring | Centring vectors |
      |---|---|---|
      | P | primitive | :math:`\\mathbf{0}` |
      | A | A-face centred | :math:`\\frac{\\mathbf{b}_s+\\mathbf{c}_s}{2}` |
      | B | B-face centred | :math:`\\frac{\\mathbf{c}_s+\\mathbf{a}_s}{2}` |
      | C | C-face centred | :math:`\\frac{\\mathbf{a}_s+\\mathbf{b}_s}{2}` |
      | I | body centred (*Innenzentriert*) | :math:`\\frac{\\mathbf{a}_s+\\mathbf{b}_s+\\mathbf{c}_s}{2}` |
      | F | all-face centred | :math:`\\frac{\\mathbf{b}_s+\\mathbf{c}_s}{2}`, :math:`\\frac{\\mathbf{c}_s+\\mathbf{a}_s}{2}`, :math:`\\frac{\\mathbf{a}_s+\\mathbf{b}_s}{2}` |
      | R | rhombohedrally centred (hexagonal axes) | :math:`\\frac{2\\mathbf{a}_s+\\mathbf{b}_s+\\mathbf{c}_s}{3}` :math:`\\frac{\\mathbf{a}_s+2\\mathbf{b}_s+2\\mathbf{c}_s}{3}` |
    
      For further details, see the `IUCr Online Dictionary of Crystallography`__.
    
      .. _website: http://reference.iucr.org/dictionary/Centred_lattice
      __ website_
      
    
    Members:
    
      invalid
    
      P : primitive
    
      A : A-face centred
    
      B : B-face centred
    
      C : C-face centred
    
      I : body-centred
    
      F : face centred
    
      R : rhombohedrally centred
    """
    A: typing.ClassVar[Bravais]
    B: typing.ClassVar[Bravais]
    C: typing.ClassVar[Bravais]
    F: typing.ClassVar[Bravais]
    I: typing.ClassVar[Bravais]
    P: typing.ClassVar[Bravais]
    R: typing.ClassVar[Bravais]
    __members__: typing.ClassVar[dict[str, Bravais]]
    invalid: typing.ClassVar[Bravais]
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class PrimitiveTransform:
    @staticmethod
    def __init__(*args, **kwargs) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def P(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def Pt(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def does_anything(self) -> bool:
        ...
    @property
    def invP(self) -> numpy.ndarray[numpy.int32]:
        ...
    @property
    def invPt(self) -> numpy.ndarray[numpy.int32]:
        ...
    @property
    def is_primitive(self) -> bool:
        ...
class Spacegroup:
    """
    The space group information as in :py:mod:`spglib:
    
      Equivalent to the struct `SpacegroupType` from
      `spg_database.h <https://github.com/spglib/spglib/blob/develop/src/spg_database.h>`_
    """
    @staticmethod
    def __init__(*args, **kwargs) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def choice(self) -> str:
        ...
    @property
    def hall_number(self) -> int:
        ...
    @property
    def hall_symbol(self) -> str:
        ...
    @property
    def international_table_full(self) -> str:
        ...
    @property
    def international_table_number(self) -> int:
        ...
    @property
    def international_table_short(self) -> str:
        ...
    @property
    def international_table_symbol(self) -> str:
        ...
    @property
    def pointgroup_number(self) -> int:
        ...
    @property
    def schoenflies_symbol(self) -> str:
        ...
class Pointgroup:
    """
    Point group information as in :py:mod:`spglib`
    
      Wrapped access to the struct originally in
      `pointgroup.h <https://github.com/spglib/spglib/blob/develop/src/pointgroup.h>`_
    """
    @staticmethod
    def __init__(*args, **kwargs) -> None:
        """
        Initialize from the serial index in the static `pointgroup_data` array
        """
    def __repr__(self) -> str:
        ...
    @property
    def holohedry(self) -> str:
        """
        Return a string representation of the Holohedry value
        
          One of `triclinic`, `monoclinic`, `orthogonal`, `tetragonal`, `trigonal`, `hexagonal`, or `cubic`.
        """
    @property
    def laue(self) -> str:
        """
        Return a string representation of the Laue class
        
          One of `1`, `2m`, `mmm`, `4m`, `4mmm`, `3`, `3m`, `6m`, `6mmm`, `m3`, or `m3m`.
        """
    @property
    def number(self) -> int:
        ...
    @property
    def symbol(self) -> str:
        ...
class SortingStatus:
    """
        An object representing the status of a single object under sorting.
    
        Internally may be represented as a single unsigned integer where two or
        more bits are reserved for various flags, or as a set of boolean values
        and an integer.
    
      
    """
    def __init__(self, sorted: bool, locked: bool, visits: int) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def locked(self) -> bool:
        """
        Return the locked flag
        """
    @property
    def sorted(self) -> bool:
        """
        Return the sorted flag
        """
    @property
    def visits(self) -> int:
        """
        Return the visit count
        """
class Symmetry:
    """
      One or more symmetry operations of a space group.
    
      A symmetry operation is the combination of a generalised rotation, :math:`W`,
      and translation, :math:`\\mathbf{w}`.
      For any position in space, :math:`\\mathbf{x}`, the operation transforms 
      :math:`\\mathbf{x}` to another equivalent position
    
      .. math::
        \\mathbf{x}' = W \\mathbf{x} + \\mathbf{w}
    
      and can equivalently be expressed as :math:`\\mathbf{x}' = \\mathscr{M}\\mathbf{x}`.
    
      Crystallographic symmetry operations have an order, :math:`o`,
      for which :math:`\\mathscr{M}^o = \\mathscr{E}` i.e. :math:`o` repeated
      applications of the operation is equivalent to the identity operator.
    
      A set of symmetry operations can form a group, :math:`\\mathbb{G}`, with the property that
      :math:`\\mathscr{M}_k = \\mathscr{M}_i \\mathscr{M}_j` with :math:`\\mathscr{M}_i,\\mathscr{M}_j,\\mathscr{M}_k \\in \\mathbb{G}`.
    
      This class can be used to hold any number of related symmetry operators, and to generate all spacegroup operators from those stored.
    
      Parameters
      ----------
      hall : int
          The integer Hall number for the desired space group operations [[deprecated]].
      W : arraylike, int
          The generalised rotation (matrix) part of the symmetry operator(s)
      w : arraylike, float
          The translation (vector) part of the symmetry operator(s)
      cifxyz : str
          The symmetry operator(s) encoded in CIF xyz format
    
      Note
      ----
      The overloaded forms of ``__init__`` take one of **hall**, (**W**, **w**), *or* **cifxyz**.
      
    """
    __hash__: typing.ClassVar[None] = None
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None:
        ...
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None:
        ...
    def __eq__(self, arg0: Symmetry) -> bool:
        ...
    @typing.overload
    def __init__(self, W: numpy.ndarray[numpy.int32], w: numpy.ndarray[numpy.float64]) -> None:
        ...
    def __len__(self) -> int:
        ...
    def generate(self) -> Symmetry:
        ...
    def generators(self) -> Symmetry:
        ...
    @property
    def W(self) -> numpy.ndarray[numpy.int32]:
        ...
    @property
    def centring(self) -> Bravais:
        ...
    @property
    def size(self) -> int:
        ...
    @property
    def w(self) -> numpy.ndarray[numpy.float64]:
        ...
class PointSymmetry:
    """
      Holds the :math:`3 \\times 3` rotation matrices :math:`R` which comprise
      a point group symmetry.
    
      A point group describes the local symmetry of a lattice point. It contains
      all of the generalised rotations of a :py:class:`~brille._brille.Symmetry`
      with none of its translations.
      
    """
    @typing.overload
    def __init__(self, Hall_number: int, time_reversal: int = 0) -> None:
        """
            Deprecated: use ``PointSymmetry(Symmetry(...))`` or :py:attr:`Lattice.pointgroup`.
        """
    @typing.overload
    def __init__(self, Symmetry: Symmetry) -> None:
        ...
    def nfolds(self, arg0: int) -> PointSymmetry:
        ...
    @property
    def W(self) -> numpy.ndarray[numpy.int32]:
        ...
    @property
    def axis(self) -> numpy.ndarray[numpy.int32]:
        ...
    @property
    def generate(self) -> PointSymmetry:
        ...
    @property
    def generators(self) -> PointSymmetry:
        ...
    @property
    def isometry(self) -> numpy.ndarray[numpy.int32]:
        ...
    @property
    def order(self) -> numpy.ndarray[numpy.int32]:
        ...
    @property
    def size(self) -> int:
        ...
class HallSymbol:
    """
        A crystallographic spacegroup's symmetries encoded in Hall's notation
    
        Hall proposed a compact unambiguous notation for the representation of the
        generators of a spacegroup. Within his notation each motion is comprised of
        a character with one or more subscripts and superscripts which describe its
        order, unique axis, and translation. The notation specifies that, depending
        on the position of a motion and details of any preceding motion, some or
        all of the sub- and superscripts can be omitted. The :class:`HallSymbol` has
        been written to handle the logic necessary to decode a Hall symbol into its
        equivalent motions.
        An added complication arises when the Hall symbol is encoded as an ASCII
        string. Namely, there are no sub- or superscript glyphs and some scheme must
        be enacted to represent them.
    
      
    """
    def __init__(self, Hall_symbol: str) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def generators(self) -> Symmetry:
        ...
class Basis:
    """
    An atom basis in a unit cell
    
      The positions and types of all symmetry-distinct atoms in a lattice define the atom basis.
      Two equivalent-type atoms may exchange position within the unit cell under application of
      a symmetry of the spacegroup.
    """
    @typing.overload
    def __init__(self, positions: numpy.ndarray[numpy.float64]) -> None:
        """
        Given only atom positions, assume all are unique types
        """
    @typing.overload
    def __init__(self, positions: numpy.ndarray[numpy.float64], types: list[int]) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def positions(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def size(self) -> int:
        ...
    @property
    def types(self) -> list[int]:
        ...
class Lattice:
    """
      A space-spanning lattice in three dimensions
    
      A space-spanning lattice in :math:`N` dimensions has :math:`N` basis vectors
      which can be described fully by their :math:`N` lengths and the
      :math:`\\sum_1^{N-1} 1` angles between each set of basis vectors, or
      :math:`\\sum_1^N 1 = \\frac{1}{2}N(N+1)` scalars in total.
      This class stores the basis vectors of the lattice described in an orthonormal space,
      plus the metric of the space, and the equivalent information for the dual of the lattice.
    
      Attributes
      ----------
      a,b,c : float
            The basis vector lengths
      alpha,beta,gamma : float
            The angles between the basis vectors, internally always in radian
      volume : float
            The volume of the lattice unit cell in units of length cubed.
      bravais : :py:class:`~brille._brille.Bravais`
            The centring type of the lattice
      spacegroup : :py:class:`~brille._brille.Symmetry`
            The Spacegroup symmetry operations of the lattice
      pointgroup : :py:class:`~brille._brille.PointSymmetry`
            The Pointgroup symmetry operations of the lattice
      basis : :py:class:`~brille._brille.Basis`
            The positions of all atoms within the lattice unit cell
      
    """
    __hash__: typing.ClassVar[None] = None
    @staticmethod
    def __eq__(*args, **kwargs):
        ...
    @staticmethod
    def __init__(*args, **kwargs):
        ...
    @staticmethod
    def __repr__(*args, **kwargs):
        ...
    @staticmethod
    def get_contravariant_metric_tensor(*args, **kwargs):
        """
          Calculate the contravariant metric tensor of the lattice
        
          Returns
          -------
          matrix_like
                The inverse of the metric of the lattice
        
                .. math::
                    g^{ij} =
                    \\begin{pmatrix}
                    a^2 & ab\\cos\\gamma & ac\\cos\\beta \\\\
                    ab\\cos\\gamma & b^2 & bc\\cos\\alpha \\\\
                    ac\\cos\\beta & bc\\cos\\alpha & c^2
                    \\end{pmatrix}^{-1}
        
          
        """
    @staticmethod
    def get_covariant_metric_tensor(*args, **kwargs):
        """
          Calculate the covariant metric tensor of the lattice
        
          Returns
          -------
          matrix_like
                The metric of the lattice
        
                .. math::
                    g_{ij} =
                    \\begin{pmatrix}
                    a^2 & ab\\cos\\gamma & ac\\cos\\beta \\\\
                    ab\\cos\\gamma & b^2 & bc\\cos\\alpha \\\\
                    ac\\cos\\beta & bc\\cos\\alpha & c^2
                    \\end{pmatrix}
        
        
          
        """
    @staticmethod
    def metric(*args, **kwargs):
        ...
    @staticmethod
    def str(*args, **kwargs):
        ...
    @staticmethod
    def vector(*args, **kwargs):
        ...
    @staticmethod
    def vectors(*args, **kwargs):
        ...
    @property
    def a(*args, **kwargs):
        ...
    @property
    def a_star(*args, **kwargs):
        ...
    @property
    def alpha(*args, **kwargs):
        ...
    @property
    def alpha_star(*args, **kwargs):
        ...
    @property
    def b(*args, **kwargs):
        ...
    @property
    def b_star(*args, **kwargs):
        ...
    @property
    def basis(*args, **kwargs):
        ...
    @property
    def beta(*args, **kwargs):
        ...
    @property
    def beta_star(*args, **kwargs):
        ...
    @property
    def bravais(*args, **kwargs):
        ...
    @property
    def c(*args, **kwargs):
        ...
    @property
    def c_star(*args, **kwargs):
        ...
    @property
    def centring_vectors(*args, **kwargs):
        """
        The centring vectors of the cell, in its fractional coordinates: the zero vector, and one
        more for each extra lattice point in a centred cell (1 for P, 2 for A, B, C and I, 3 for
        R, 4 for F).
        """
    @property
    def gamma(*args, **kwargs):
        ...
    @property
    def gamma_star(*args, **kwargs):
        ...
    @property
    def pointgroup(*args, **kwargs):
        ...
    @property
    def primitive_basis(*args, **kwargs):
        """
        The atoms of one primitive cell: the first given atom of each centring orbit of
        :py:attr:`basis`, at its given position.
        
        Eigenvectors given to a grid describe these atoms, in this order. For a primitive cell
        this is :py:attr:`basis` itself; for a centred conventional cell it is a half, a third
        or a quarter of it (see :py:func:`brille.utils.conventional_to_primitive`).
        """
    @property
    def real_vectors(*args, **kwargs):
        ...
    @property
    def reciprocal_vectors(*args, **kwargs):
        ...
    @property
    def spacegroup(*args, **kwargs):
        ...
    @spacegroup.setter
    def spacegroup(*args, **kwargs):
        ...
    @property
    def volume(*args, **kwargs):
        ...
    @property
    def volume_star(*args, **kwargs):
        ...
class Polyhedron:
    @typing.overload
    def __init__(self, vertices: numpy.ndarray[numpy.float64]) -> None:
        ...
    @typing.overload
    def __init__(self, vertices: numpy.ndarray[numpy.float64], faces: list[list[int]]) -> None:
        ...
    def intersection(self, arg0: Polyhedron) -> Polyhedron:
        ...
    @property
    def centre(self) -> Polyhedron:
        ...
    @property
    def faces(self) -> list[list[int]]:
        ...
    @property
    def mirror(self) -> Polyhedron:
        ...
    @property
    def normals(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def points(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def vertices(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def volume(self) -> float:
        ...
class LPolyhedron:
    def intersection(self, arg0: LPolyhedron) -> LPolyhedron:
        ...
    def rotate(self, arg0: PointSymmetry, arg1: int) -> LPolyhedron:
        ...
    def to_Cartesian(self) -> Polyhedron:
        ...
    @property
    def centre(self) -> LPolyhedron:
        ...
    @property
    def faces(self) -> list[list[int]]:
        ...
    @property
    def mirror(self) -> LPolyhedron:
        ...
    @property
    def normals(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def points(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def vertices(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def volume(self) -> float:
        ...
class BrillouinZone:
    """
        Construct and hold a first Brillouin zone and, optionally and by default,
        an irreducible Brillouin zone.
    
        The region closer to a given lattice point than to any other is the
        Wigner-Seitz cell of that lattice. The same construction is one possible
        first Brillouin zone of a reciprocal lattice and is used within ``brille``.
        For example, a two-dimensional hexagonal lattice has a first Brillouin
        zone which is a hexagon:
    
        .. tikz::
            :libs: calc
    
            \\begin{tikzpicture}[scale=5,dot/.style = {fill,radius=0.02},
            latpt/.style = {color=gray!50!white, fill},]
            \\coordinate (astar) at (1,0);
            \\coordinate (bstar) at (0.5,0.8660254037844386);
            \\clip ($-0.1*(bstar)$) rectangle ($3*(astar)+2.1*(bstar)$);
            %
            \\coordinate (C1) at ($0.3333*(astar)+0.3333*(bstar)$);
            \\coordinate (C2) at ($-0.3333*(astar)+0.6667*(bstar)$);
            \\coordinate (C3) at ($-0.6667*(astar)+0.3333*(bstar)$);
            \\coordinate (C4) at ($-0.3333*(astar)-0.3333*(bstar)$);
            \\coordinate (C5) at ($0.3333*(astar)-0.6667*(bstar)$);
            \\coordinate (C6) at ($0.6667*(astar)-0.3333*(bstar)$);
            %
            \\foreach \\h in {-1,...,6}{
            \\foreach \\k in {-1,...,3}{
              \\draw[latpt] ($\\h*(astar) + \\k*(bstar)$) circle[dot];
              \\draw[color=yellow!75!black,dotted] ($\\h*(astar)+\\k*(bstar)$) +(C1)
                -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
            }}
            %
            \\coordinate (G) at ($(astar)+(bstar)$);
            \\draw[color=yellow!75!black, line width=1mm] (G) +(C1) -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
            \\end{tikzpicture}
    
    
        Since all physical properties of a crystal must have the same periodicity
        as its lattice, the powerful feature of the first Brillouin zone is that it
        encompasses a region of reciprocal space which must fully represent all
        of reciprocal space.
    
        Most crystals contain rotational or rotoinversion symmetries in addition to
        the translational ones which give rise to the first Brillouin zone. These
        symmetries are the pointgroup of the lattice and enforce that the properties
        of the crystal also have the same symmetry. The first Brillouin zone,
        therefore, typically contains redundant information.
    
        An irreducible Brillouin zone is a subsection of the first Brillouin zone
        which contains the minimal part required to have only unique crystal
        properties. This class can find an irreducible Brillouin zone for any
        crystal lattice. In the example of the hexagonal lattice there are six
        equivalent irreducible Brillouin zones one of which is:
    
        .. tikz::
            :libs: calc
    
            \\begin{tikzpicture}[scale=5,dot/.style = {fill,radius=0.02},
            latpt/.style = {color=gray!50!white, fill},]
            \\coordinate (astar) at (1,0);
            \\coordinate (bstar) at (0.5,0.8660254037844386);
            \\clip ($-0.1*(bstar)$) rectangle ($3*(astar)+2.1*(bstar)$);
            %
            \\coordinate (C1) at ($0.3333*(astar)+0.3333*(bstar)$);
            \\coordinate (C2) at ($-0.3333*(astar)+0.6667*(bstar)$);
            \\coordinate (C3) at ($-0.6667*(astar)+0.3333*(bstar)$);
            \\coordinate (C4) at ($-0.3333*(astar)-0.3333*(bstar)$);
            \\coordinate (C5) at ($0.3333*(astar)-0.6667*(bstar)$);
            \\coordinate (C6) at ($0.6667*(astar)-0.3333*(bstar)$);
            %
            \\foreach \\h in {-1,...,6}{
            \\foreach \\k in {-1,...,3}{
              \\draw[latpt] ($\\h*(astar) + \\k*(bstar)$) circle[dot];
              \\draw[color=yellow!75!black,dotted] ($\\h*(astar)+\\k*(bstar)$) +(C1)
                -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
            }}
            %
            \\coordinate (G) at ($(astar)+(bstar)$);
            \\draw[color=yellow!75!black,line width=1mm] (G) -- +(C4) -- +(C5) -- cycle;
            \\end{tikzpicture}
    
        Parameters
        ----------
        lattice: :py:class:`brille._brille.Reciprocal`
            The reciprocal space lattice for which a Brillouin zone will be found
        use_primitive: bool
            If the provided :py:class:`brille._brille.Reciprocal` lattice is a
            conventional Bravais lattice, this parameter controls whether the
            equivalent primitive Bravais lattice should be used to find the first
            Brillouin zone. This is ``True`` by default and should only be modified
            for testing purposes.
        search_length: int
            The Wigner-Seitz construction of the first Brillouin zone finds the
            volume of space closer to a chosen reciprocal lattice point than any
            other reciprocal lattice point. This is accomplished by successively
            dividing the space by planes halfway between the chosen point and a
            subset of all other planes. The subset used is controlled by
            `search_length` and is every unique :math:`(\\pm s_i\\,0\\,0)`,
            :math:`(0\\,\\pm s_j\\,0)`, :math:`(0\\,0\\,\\pm s_k)`,
            :math:`(\\pm s_i\\,\\pm s_j\\,0)`, :math:`(\\pm s_i\\,0\\,\\pm s_k)`,
            :math:`(0\\,\\pm s_j\\,\\pm s_k)`, :math:`(\\pm s_i\\,\\pm s_j\\,\\pm s_k)` for
            :math:`1 \\le s_\\alpha \\le` `search_length`.
            If the reciprocal lattice is primitive then the default `search_length`
            of ``1`` should always give the correct first Brillouin zone.
            For extra assurance that the correct first Brillouin zone is found, the
            procedure is internally repeated with `search_length` incremented by
            one and an error is raised if the two constructed polyhedra have
            different volumes.
        time_reversal_symmetry: bool
            Controls whether time reversal symmetry should be added to pointgroups
            lacking space inversion. This affects the found irreducible Brillouin
            zone for such systems. To avoid inadvertently adding time reversal
            symmetry when it is not appropriate, this is ``False`` by default.
            Time reversal holds for non-magnetic systems. It is anti-unitary: phonon
            eigenvectors (``RotatesLike.Gamma``) at a time-reversed point
            :math:`-R\\mathbf{q}` are the complex conjugates of those at
            :math:`R\\mathbf{q}`, so no inversion symmetry of the crystal is needed.
        wedge_search: bool
            Controls whether an irreducible Brillouin zone should be found. With
            this set to ``False`` the returned :py:class:`brille._brille.BrillouinZone`
            will only contain the first Brillouin zone. If ``True`` the pointgroup
            symmetry operations will be used to identify *an* irreducible Brillouin
            zone as well. If the provided lattice's parameters do not match the
            symmetry of the pointgroup (e.g., a lattice which should be tetragonal
            like :math:`I4/mmm` but constructed with :math:`\\gamma=120^\\circ`) the
            algorithm will fail to find an appropriate irreducible Brillouin zone
            and an error will be raised. (Set to ``True`` by default).
        warn_near_symmetry: bool
            Keyword only. Whether to warn, with a
            :py:class:`brille.NearSymmetryWarning`, when the lattice is within
            1e-4 of a lattice with more symmetry operations (e.g., rhombohedral
            but nearly cubic). Such a zone has faces or edges much smaller than
            itself, which make meshes slow and poorly shaped. The zone is not
            changed. (Set to ``True`` by default).
      
    """
    @staticmethod
    def from_file(filename: str, entry: str = 'BrillouinZone') -> BrillouinZone:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/bz", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    @typing.overload
    def __init__(self, lattice: Lattice, use_primitive: bool = True, search_length: int = 1, time_reversal_symmetry: bool = False, wedge_search: bool = True, divide_primitive: bool = True, *, warn_near_symmetry: bool = True) -> None:
        ...
    @typing.overload
    def __init__(self, lattice: Lattice, approx_config: ApproxConfig, use_primitive: bool = True, search_length: int = 1, time_reversal_symmetry: bool = False, wedge_search: bool = True, divide_primitive: bool = True, *, warn_near_symmetry: bool = True) -> None:
        ...
    def ir_moveinto(self, Q: numpy.ndarray[numpy.float64], threads: int = 0) -> tuple:
        """
            Find points equivalent to those provided within the irreducible Brillouin zone.
        
            The BrillouinZone object defines a volume of reciprocal space which contains
            an irreducible part of the full reciprocal-space. This method will find
            points equivalent under the operations of the lattice which fall within this
            irreducible volume.
        
            Parameters
            ----------
            Q : :py:class:`numpy.ndarray`
                A 2 dimensional array of three-vectors (``Q.shape[1]==3``) expressed in
                units of the reciprocal lattice.
            threads : integer, optional
                The number of parallel threads that should be used. If this value is less
                than one, the ``BRILLE_NUM_THREADS`` environment variable sets the number,
                or one thread per logical core is used if it is not set.
        
            Returns
            -------
            Qir : :py:class:`numpy.ndarray`
                The array of equivalent irreducible :math:`\\mathbf{q}_\\text{ir}` points
                for all :math:`\\mathbf{Q}`;
            tau : :py:class:`numpy.ndarray`
                the closest reciprocal lattice vector, :math:`\\boldsymbol{\\tau}`,
                to each :math:`\\mathbf{Q}`;
            R : :py:class:`numpy.ndarray`
                the pointgroup symmetry operation :math:`R`
            Rinv : :py:class:`numpy.ndarray`
                the inverse point group symmetry operation which obey
                :math:`\\mathbf{Q} = R^{-1} \\mathbf{q}_\\text{ir} + \\boldsymbol{\\tau}`.
        """
    def ir_moveinto_wedge(self, Q: numpy.ndarray[numpy.float64], threads: int = 0) -> tuple:
        """
            Find points equivalent to those provided within the irreducible wedge.
        
            The BrillouinZone object defines a wedge of reciprocal space which contains
            an irreducible part of the full-space 4π steradian solid angle. This method
            will find points equivalent under the pointgroup operations of the lattice
            which fall within this irreducible solid angle and maintain their absolute
            magnitude.
        
            Parameters
            ----------
            Q : :py:class:`numpy.ndarray`
                A 2 dimensional array of three-vectors (``Q.shape[1]==3``) expressed in
                units of the reciprocal lattice.
            threads : integer, optional (default 0)
                The number of parallel threads that should be used. If this value is less
                than one, the ``BRILLE_NUM_THREADS`` environment variable sets the number,
                or one thread per logical core is used if it is not set.
        
            Returns
            -------
            :py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`
                The array of equivalent in-wedge :math:`\\mathbf{Q}_\\text{ir}` points
                for all :math:`\\mathbf{Q}`, and the pointgroup operation fulfilling
                :math:`\\mathbf{Q}_\\text{ir} = R \\mathbf{Q}`.
        """
    def isinside(self, points: numpy.ndarray[numpy.float64]) -> list[bool]:
        """
            Determine whether each of the provided reciprocal lattice points is located
            within the first Brillouin zone
        
            Parameters
            ----------
            Q : :py:class:`numpy.ndarray`
                A 2 dimensional array of three-vectors (``Q.shape[1]==3``) expressed in
                units of the reciprocal lattice.
        
            Returns
            -------
            :py:class:`numpy.ndarray`
                One dimensional logical array with ``True`` indicating 'inside'
        """
    def lattice_symmetry_counts(self, exact: float = 1e-10, near: float = 0.0001) -> tuple[int, int]:
        """
          The number of symmetry operations of the lattice itself, exactly and nearly
        
          Counts the lattice's own symmetries (its holohedry, not the crystal's
          symmetry) within `exact` and within `near`, relative to the lattice metric.
        
          Returns
          -------
          tuple[int, int]
              The counts within `exact` and within `near`; more near than exact
              symmetries mean the lattice is close to a more symmetric one.
        """
    def moveinto(self, Q: numpy.ndarray[numpy.float64], threads: int = 0) -> tuple:
        """
            Find points equivalent to those provided within the first Brillouin zone.
        
            Parameters
            ----------
            Q : :py:class:`numpy.ndarray`
                A 2 dimensional array of three-vectors (``Q.shape[1]==3``) expressed in
                units of the reciprocal lattice.
            threads : integer, optional
                The number of parallel threads that should be used. If this value is less
                than one, the ``BRILLE_NUM_THREADS`` environment variable sets the number,
                or one thread per logical core is used if it is not set.
        
            Returns
            -------
            :py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`
                The floating point array of equivalent reduced :math:`\\mathbf{q}`
                points for all :math:`\\mathbf{Q}`, and an integer array filled with
                :math:`\\boldsymbol{\\tau} = \\mathbf{Q}-\\mathbf{q}`.
        """
    def to_file(self, filename: str, entry: str = 'BrillouinZone', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/bz", where to write inside the file,
              with a default equal to BrillouinZone name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def faces_per_vertex(self) -> list[list[int]]:
        """
          Return the first Brillouin zone face indices for each unique face corner
        """
    @property
    def half_edge_points(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone face edge centres in rlu
        """
    @property
    def half_edge_points_invA(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone face edge centres in inverse ångstrom
        """
    @property
    def ir_faces_per_vertex(self) -> list[list[int]]:
        """
          Return the irreducible Brillouin zone face index per unique face corner
        """
    @property
    def ir_normals(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the irreducible Brillouin zone face normals in rlu
        """
    @property
    def ir_normals_invA(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the irreducible Brillouin zone face normals in inverse ångstrom
        """
    @property
    def ir_normals_primitive(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the irreducible Brillouin zone face normals in primitive-lattice rlu
        """
    @property
    def ir_points(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the irreducible Brillouin zone face centres in rlu
        """
    @property
    def ir_points_invA(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the irreducible Brillouin zone face centres in inverse ångstrom
        """
    @property
    def ir_points_primitive(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the irreducible Brillouin zone face centres in primitive-lattice rlu
        """
    @property
    def ir_polyhedron(self) -> LPolyhedron:
        """
          Returns the irreducible Brillouin zone :py:class:`brille._brille.Polyhedron`
        
          Returns
          -------
          :py:class:`brille._brille.Polyhedron`
              If no irreducible Brillouin zone was requested at construction, the returned
              polyhedron is that of the first Brillouin zone instead.
        """
    @property
    def ir_polyhedron_generated(self) -> LPolyhedron:
        """
          Returns the found irreducible Brillouin zone :py:class:`brille._brille.Polyhedron`
        
          If the lattice pointgroup does not contain the space inversion operator
          the internally held 'irreducible' polyhedron is only half of the real
          irreducible polyhedron. This method gives access to the polyhedron found by
          the algorithm before being doubled for output.
        """
    @property
    def ir_vertices(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the irreducible Brillouin zone unique face corners in rlu
        """
    @property
    def ir_vertices_invA(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the irreducible Brillouin zone unique face corners in inverse ångstrom
        """
    @property
    def ir_vertices_per_face(self) -> list[list[int]]:
        """
          Return the irreducible Brillouin zone unique face corners per face
        """
    @property
    def ir_vertices_primitive(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the irreducible Brillouin zone unique face corners in primitive-lattice rlu
        """
    @property
    def lattice(self) -> Lattice:
        """
          Returns the defining :py:class:`brille._brille.Lattice` lattice
        """
    @property
    def normals(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone face normals in rlu
        """
    @property
    def normals_invA(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone face normals in inverse ångstrom
        """
    @property
    def normals_primitive(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone face normals in primitive-lattice rlu
        """
    @property
    def points(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone face centres in rlu
        """
    @property
    def points_invA(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone face centres in inverse ångstrom
        """
    @property
    def points_primitive(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone face centres in primitive-lattice rlu
        """
    @property
    def polyhedron(self) -> LPolyhedron:
        """
          Returns the first Brillouin zone :py:class:`brille._brille.Polyhedron`
        """
    @property
    def vertices(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone unique face corners in rlu
        """
    @property
    def vertices_invA(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone unique face corners in inverse ångstrom
        """
    @property
    def vertices_per_face(self) -> list[list[int]]:
        """
          Return the first Brillouin zone face corner indices for each face
        """
    @property
    def vertices_primitive(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the first Brillouin zone unique face corners in primitive-lattice rlu
        """
    @property
    def wedge_normals(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the normals of the irreducible wedge rlu
        """
    @property
    def wedge_normals_invA(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the normals of the irreducible wedge inverse ångstrom
        """
    @property
    def wedge_normals_primitive(self) -> numpy.ndarray[numpy.float64]:
        """
          Return the normals of the irreducible wedge primitive-lattice rlu
        """
class BZMeshQdd:
    @staticmethod
    def from_file(filename: str, entry: str = 'BZMeshQdd') -> BZMeshQdd:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/grid", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    def __init__(self, brillouin_zone: BrillouinZone, max_size: float = -1.0, num_levels: int = 3, max_points: int = -1) -> None:
        """
        A structured tetrahedral mesh of a Brillouin zone's irreducible part
        
        A grid of the reciprocal lattice, divided finely enough for ``max_size``, is
        clipped exactly to the irreducible zone.
        
        Parameters
        ----------
        brillouin_zone : BrillouinZone
            The zone whose irreducible part the mesh fills.
        max_size : float, optional (default: -1)
            The largest tetrahedron volume, in cubic reciprocal Angstrom, which sets the
            grid spacing; if not positive, the grid is the reciprocal lattice itself.
            Each grid cell holds six tetrahedra, so ``max_size = node_volume_fraction /
            6`` gives about as many vertices as a :py:class:`BZTrellisQdc` with that
            ``node_volume_fraction``, and ``brillouin_zone.ir_polyhedron.volume / (6 *
            points)`` gives roughly 1.5 to 3 times ``points`` vertices, the most for
            small meshes.
        num_levels : int, optional
            Unused; kept for compatibility.
        max_points : int, optional (default: -1)
            If positive, the grid is coarsened until its estimated number of vertices is
            at most this, with a RuntimeWarning (see :py:attr:`refinement_limited`).
        """
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], vectors_data: numpy.ndarray[numpy.float64], vectors_elements: numpy.ndarray[numpy.int32], sort: bool = False) -> None:
        """
        Provide data required for interpolation to the grid without cost information.
        
        .. Note
        .. ----
        .. This method should probably be followed by :py:meth:`set_cost_info` prior to
        .. any attempt to interpolate the data in the grid.
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below).
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below).
        
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            ``values_data``
        vectors_elements:
            Like ``values_elements`` but for the eigenvectors
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_data: numpy.ndarray[numpy.float64], vectors_elements: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
        Provide all data required for interpolation to the grid at once
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below),
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below)
            * which scalar cost function should be used (see below),
            * which vector cost function should be used (see below).
        
            See the note below for the meaning of the last three values.
        values_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvalue elements stored in the grid
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            **values_data**
        vectors_elements:
            Like **values_elements** but for the eigenvectors
        vectors_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvector elements stored in the grid
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        
          Mapping of integers to scalar cost function:
        
          | value | function(x,y) |
          |---|---|
          | 0 | magnitude(x-y) |
        
          Mapping of integers to vector cost function:
        
          | value | function(vec_x, vec_y) |
          |---|---|
          | 0 | sin(hermitian_angle(vec_x, vec_y)) |
          | 1 | vector_distance(vec_x, vec_y) |
          | 2 | 1 - vector_product(vec_x, vec_y) |
          | 3 | vector_angle(vec_x, vec_y) |
          | 4 | hermitian_angle(vec_x, vec_y) |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def ir_interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.float64], numpy.ndarray[numpy.float64]]:
        """
          Perform linear interpolation of the stored data at irreducible equivalent points
        
          The irreducible first Brillouin zone is the part of reciprocal space which is
          invariant under application of the integer translations *and* the pointgroup
          operations of a reciprocal space lattice. This method finds points equivalent
          to the input within the irreducible first Brillouin zone and then interpolates
          pre-stored information to provide an estimate at the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              irreducible first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def refine(self, where: typing.Any = None, values: typing.Any = None, vectors: typing.Any = None, resolution: float | None = None, points_per_resolution: float = 2.0) -> numpy.ndarray[numpy.float64]:
        """
        Refine the mesh by bisecting tetrahedra; existing vertices keep their indices and data.
        
        Parameters
        ----------
        where, resolution, points_per_resolution
            As for :py:meth:`refinement_points`, which gives the points this adds.
        values, vectors : numpy.ndarray, optional
            If the mesh holds data (after :py:meth:`fill`), the data for the new vertices,
            laid out per point as the filled data and for exactly the points
            :py:meth:`refinement_points` returns, in that order. Not allowed before the
            mesh is filled.
        
        Returns
        -------
        numpy.ndarray
            The new vertices, shape (N, 3), in relative lattice units, appended to
            :py:attr:`rlu` in this order.
        
        Note
        ----
        The mode permutations found by :py:meth:`sort` are reset; sort again after
        refining if needed.
        """
    def refinement_points(self, where: typing.Any = None, resolution: float | None = None, points_per_resolution: float = 2.0) -> numpy.ndarray[numpy.float64]:
        """
        The points that :py:meth:`refine` would add, without changing the mesh.
        
        Evaluate your model at these points and compare with
        :py:meth:`ir_interpolate_at` there to decide whether refining is worth it; then
        pass the model's values to :py:meth:`refine` with the same arguments.
        
        Parameters
        ----------
        where : None, bool array or int array
            The tetrahedra to split: all of them (None), those where a boolean mask with
            one entry per tetrahedron is true, or those with the given indices.
            Neighbouring tetrahedra are split as needed to keep the mesh conforming, and
            split edges on the zone boundary are split with their symmetry equivalents, so
            that equivalent zone faces keep matching.
        resolution : float, optional
            The resolution limit, in inverse Angstrom. No edge is split to below
            ``resolution / points_per_resolution``: a tetrahedron whose longest edge is at
            most twice that is left whole.
        points_per_resolution : float, optional (default: 2)
            How finely to resolve ``resolution``.
        
        Returns
        -------
        numpy.ndarray
            The new vertices, shape (N, 3), in relative lattice units like :py:attr:`rlu`.
            :py:meth:`refine` appends them to the vertices in this order.
        """
    def release_triangulation(self) -> None:
        """
        Free the memory refinement holds between refinements.
        
        After :py:meth:`refine` (or :py:meth:`refinement_points`) the mesh keeps the
        triangulation refinement works on, several times the memory of the mesh itself.
        This frees it. The mesh is unchanged and can still be refined: the triangulation
        is rebuilt when next needed, which costs a build of the mesh plus a replay of the
        refinements made so far.
        """
    def set_flags_weights(self, values_flags: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_flags: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
          Set :py:class:`~brille._brille.RotatesLike`, :py:class:`~brille._brille.LengthUnit`
          and cost functions plus relative cost weights for the values and vectors
          stored in the object
        
          Parameters
          ----------
          values_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvalues at
              neighbouring grid points for scalar- and vector-like eigenvalues.
          values_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvalue elements stored in the grid
          vectors_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvectors at
              neighbouring grid points for scalar- and vector-like eigenvectors.
          vectors_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvector elements stored in the grid
          sort : bool, optional
              Whether the equivalent-mode permutations should be (re)determined following
              the update to the flags and weights.
        
        
          Note
          ----
            Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:
        
            | value | :py:class:`RotatesLike` |
            |---|---|
            | 0 | `vector` |
            | 1 | `pseudovector` |
            | 2 | `Gamma` |
          
            Mapping of integers to :py:class:`LengthUnit` values:
        
            | value | :py:class:`LengthUnit` |
            |---|---|
            | 0 | `none` |
            | 1 | `angstrom` |
            | 2 | `inverse_angstrom` |
            | 3 | `real_lattice` |
            | 4 | `reciprocal_lattice` |
        
            Mapping of integers to scalar cost function:
          
            | value | function(x,y) |
            |---|---|
            | 0 | magnitude(x-y) |
          
            Mapping of integers to vector cost function:
          
            | value | function(vec_x, vec_y) |
            |---|---|
            | 0 | sin(hermitian_angle(vec_x, vec_y)) |
            | 1 | vector_distance(vec_x, vec_y) |
            | 2 | 1 - vector_product(vec_x, vec_y) |
            | 3 | vector_angle(vec_x, vec_y) |
            | 4 | hermitian_angle(vec_x, vec_y) |
        
            Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def set_vector_normalization(self, normalize: bool | None = True, metric: list[float] | None = None) -> None:
        """
            Choose when interpolated eigenvectors are scaled to unit norm
        
            Linear interpolation between unit eigenvectors gives vectors shorter than
            one wherever neighbouring eigenvectors differ, so structure factors computed
            from them come out too small. Normalization scales each interpolated branch
            :math:`v` to :math:`v/\\sqrt{|\\langle v|M|v\\rangle|}`.
        
            By default it is automatic: eigenvectors stored in Cartesian units
            (:py:class:`LengthUnit` ``angstrom`` or ``inverse_angstrom``, as Euphonic
            stores them) are normalized, and those in lattice units, whose length depends
            on the lattice, are not. The choice survives :py:meth:`fill` and saving to HDF5.
        
            Parameters
            ----------
            normalize : bool or None, optional
                ``True`` always normalizes, and raises a RuntimeError for eigenvectors in
                lattice units; ``False`` never does; ``None`` restores the automatic default.
            metric : float, vector-like, optional
                A diagonal metric :math:`M`, one weight per element of a branch (for
                phonons, :math:`3N`). The default is the identity, the ordinary norm.
                For Bogoliubov (spin-wave) vectors use
                :math:`\\eta=\\mathrm{diag}(1,\\ldots,1,-1,\\ldots,-1)`; the sign of
                :math:`\\langle v|\\eta|v\\rangle` is kept.
        """
    def sort(self) -> None:
        ...
    def to_file(self, filename: str, entry: str = 'BZMeshQdd', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/grid", where to write inside the file,
              with a default equal to the object Class name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def BrillouinZone(self) -> BrillouinZone:
        ...
    @property
    def bytes_per_point(self) -> int:
        """
            Return the memory required per interpolation point *result* in bytes
        """
    @property
    def holds_triangulation(self) -> bool:
        """
        Whether the triangulation that refinement works on is in memory (see release_triangulation)
        """
    @property
    def invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def normalizes_vectors(self) -> bool:
        """
            Whether interpolated eigenvectors are normalized, given the stored data; see :py:meth:`set_vector_normalization`
        """
    @property
    def refinable(self) -> bool:
        """
        Whether the mesh can be refined; a mesh read from a file written before refinement existed can't be
        """
    @property
    def refinement_limited(self) -> bool:
        """
        Whether max_points made the mesh coarser than max_size asked for
        """
    @property
    def rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def tetrahedra(self) -> numpy.ndarray[numpy.uint32]:
        ...
    @property
    def values(self) -> numpy.ndarray[numpy.float64]:
        """
            Return a shared view of the stored eigenvalues
        """
    @property
    def vector_metric(self) -> list[float]:
        """
            The diagonal metric used to normalize eigenvectors; empty for the identity
        """
    @property
    def vector_normalization(self) -> str:
        """
            When eigenvectors are normalized: ``"automatic"`` (the default), ``"on"`` or ``"off"``
        """
    @property
    def vectors(self) -> numpy.ndarray[numpy.float64]:
        """
            Return a shared view of the stored eigenvectors
        """
class BZMeshQdc:
    @staticmethod
    def from_file(filename: str, entry: str = 'BZMeshQdc') -> BZMeshQdc:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/grid", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    def __init__(self, brillouin_zone: BrillouinZone, max_size: float = -1.0, num_levels: int = 3, max_points: int = -1) -> None:
        """
        A structured tetrahedral mesh of a Brillouin zone's irreducible part
        
        A grid of the reciprocal lattice, divided finely enough for ``max_size``, is
        clipped exactly to the irreducible zone.
        
        Parameters
        ----------
        brillouin_zone : BrillouinZone
            The zone whose irreducible part the mesh fills.
        max_size : float, optional (default: -1)
            The largest tetrahedron volume, in cubic reciprocal Angstrom, which sets the
            grid spacing; if not positive, the grid is the reciprocal lattice itself.
            Each grid cell holds six tetrahedra, so ``max_size = node_volume_fraction /
            6`` gives about as many vertices as a :py:class:`BZTrellisQdc` with that
            ``node_volume_fraction``, and ``brillouin_zone.ir_polyhedron.volume / (6 *
            points)`` gives roughly 1.5 to 3 times ``points`` vertices, the most for
            small meshes.
        num_levels : int, optional
            Unused; kept for compatibility.
        max_points : int, optional (default: -1)
            If positive, the grid is coarsened until its estimated number of vertices is
            at most this, with a RuntimeWarning (see :py:attr:`refinement_limited`).
        """
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], sort: bool = False) -> None:
        """
        Provide data required for interpolation to the grid without cost information.
        
        .. Note
        .. ----
        .. This method should probably be followed by :py:meth:`set_cost_info` prior to
        .. any attempt to interpolate the data in the grid.
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below).
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below).
        
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            ``values_data``
        vectors_elements:
            Like ``values_elements`` but for the eigenvectors
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
        Provide all data required for interpolation to the grid at once
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below),
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below)
            * which scalar cost function should be used (see below),
            * which vector cost function should be used (see below).
        
            See the note below for the meaning of the last three values.
        values_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvalue elements stored in the grid
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            **values_data**
        vectors_elements:
            Like **values_elements** but for the eigenvectors
        vectors_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvector elements stored in the grid
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        
          Mapping of integers to scalar cost function:
        
          | value | function(x,y) |
          |---|---|
          | 0 | magnitude(x-y) |
        
          Mapping of integers to vector cost function:
        
          | value | function(vec_x, vec_y) |
          |---|---|
          | 0 | sin(hermitian_angle(vec_x, vec_y)) |
          | 1 | vector_distance(vec_x, vec_y) |
          | 2 | 1 - vector_product(vec_x, vec_y) |
          | 3 | vector_angle(vec_x, vec_y) |
          | 4 | hermitian_angle(vec_x, vec_y) |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def ir_interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.float64], numpy.ndarray[numpy.complex128]]:
        """
          Perform linear interpolation of the stored data at irreducible equivalent points
        
          The irreducible first Brillouin zone is the part of reciprocal space which is
          invariant under application of the integer translations *and* the pointgroup
          operations of a reciprocal space lattice. This method finds points equivalent
          to the input within the irreducible first Brillouin zone and then interpolates
          pre-stored information to provide an estimate at the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              irreducible first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def refine(self, where: typing.Any = None, values: typing.Any = None, vectors: typing.Any = None, resolution: float | None = None, points_per_resolution: float = 2.0) -> numpy.ndarray[numpy.float64]:
        """
        Refine the mesh by bisecting tetrahedra; existing vertices keep their indices and data.
        
        Parameters
        ----------
        where, resolution, points_per_resolution
            As for :py:meth:`refinement_points`, which gives the points this adds.
        values, vectors : numpy.ndarray, optional
            If the mesh holds data (after :py:meth:`fill`), the data for the new vertices,
            laid out per point as the filled data and for exactly the points
            :py:meth:`refinement_points` returns, in that order. Not allowed before the
            mesh is filled.
        
        Returns
        -------
        numpy.ndarray
            The new vertices, shape (N, 3), in relative lattice units, appended to
            :py:attr:`rlu` in this order.
        
        Note
        ----
        The mode permutations found by :py:meth:`sort` are reset; sort again after
        refining if needed.
        """
    def refinement_points(self, where: typing.Any = None, resolution: float | None = None, points_per_resolution: float = 2.0) -> numpy.ndarray[numpy.float64]:
        """
        The points that :py:meth:`refine` would add, without changing the mesh.
        
        Evaluate your model at these points and compare with
        :py:meth:`ir_interpolate_at` there to decide whether refining is worth it; then
        pass the model's values to :py:meth:`refine` with the same arguments.
        
        Parameters
        ----------
        where : None, bool array or int array
            The tetrahedra to split: all of them (None), those where a boolean mask with
            one entry per tetrahedron is true, or those with the given indices.
            Neighbouring tetrahedra are split as needed to keep the mesh conforming, and
            split edges on the zone boundary are split with their symmetry equivalents, so
            that equivalent zone faces keep matching.
        resolution : float, optional
            The resolution limit, in inverse Angstrom. No edge is split to below
            ``resolution / points_per_resolution``: a tetrahedron whose longest edge is at
            most twice that is left whole.
        points_per_resolution : float, optional (default: 2)
            How finely to resolve ``resolution``.
        
        Returns
        -------
        numpy.ndarray
            The new vertices, shape (N, 3), in relative lattice units like :py:attr:`rlu`.
            :py:meth:`refine` appends them to the vertices in this order.
        """
    def release_triangulation(self) -> None:
        """
        Free the memory refinement holds between refinements.
        
        After :py:meth:`refine` (or :py:meth:`refinement_points`) the mesh keeps the
        triangulation refinement works on, several times the memory of the mesh itself.
        This frees it. The mesh is unchanged and can still be refined: the triangulation
        is rebuilt when next needed, which costs a build of the mesh plus a replay of the
        refinements made so far.
        """
    def set_flags_weights(self, values_flags: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_flags: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
          Set :py:class:`~brille._brille.RotatesLike`, :py:class:`~brille._brille.LengthUnit`
          and cost functions plus relative cost weights for the values and vectors
          stored in the object
        
          Parameters
          ----------
          values_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvalues at
              neighbouring grid points for scalar- and vector-like eigenvalues.
          values_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvalue elements stored in the grid
          vectors_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvectors at
              neighbouring grid points for scalar- and vector-like eigenvectors.
          vectors_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvector elements stored in the grid
          sort : bool, optional
              Whether the equivalent-mode permutations should be (re)determined following
              the update to the flags and weights.
        
        
          Note
          ----
            Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:
        
            | value | :py:class:`RotatesLike` |
            |---|---|
            | 0 | `vector` |
            | 1 | `pseudovector` |
            | 2 | `Gamma` |
          
            Mapping of integers to :py:class:`LengthUnit` values:
        
            | value | :py:class:`LengthUnit` |
            |---|---|
            | 0 | `none` |
            | 1 | `angstrom` |
            | 2 | `inverse_angstrom` |
            | 3 | `real_lattice` |
            | 4 | `reciprocal_lattice` |
        
            Mapping of integers to scalar cost function:
          
            | value | function(x,y) |
            |---|---|
            | 0 | magnitude(x-y) |
          
            Mapping of integers to vector cost function:
          
            | value | function(vec_x, vec_y) |
            |---|---|
            | 0 | sin(hermitian_angle(vec_x, vec_y)) |
            | 1 | vector_distance(vec_x, vec_y) |
            | 2 | 1 - vector_product(vec_x, vec_y) |
            | 3 | vector_angle(vec_x, vec_y) |
            | 4 | hermitian_angle(vec_x, vec_y) |
        
            Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def set_vector_normalization(self, normalize: bool | None = True, metric: list[float] | None = None) -> None:
        """
            Choose when interpolated eigenvectors are scaled to unit norm
        
            Linear interpolation between unit eigenvectors gives vectors shorter than
            one wherever neighbouring eigenvectors differ, so structure factors computed
            from them come out too small. Normalization scales each interpolated branch
            :math:`v` to :math:`v/\\sqrt{|\\langle v|M|v\\rangle|}`.
        
            By default it is automatic: eigenvectors stored in Cartesian units
            (:py:class:`LengthUnit` ``angstrom`` or ``inverse_angstrom``, as Euphonic
            stores them) are normalized, and those in lattice units, whose length depends
            on the lattice, are not. The choice survives :py:meth:`fill` and saving to HDF5.
        
            Parameters
            ----------
            normalize : bool or None, optional
                ``True`` always normalizes, and raises a RuntimeError for eigenvectors in
                lattice units; ``False`` never does; ``None`` restores the automatic default.
            metric : float, vector-like, optional
                A diagonal metric :math:`M`, one weight per element of a branch (for
                phonons, :math:`3N`). The default is the identity, the ordinary norm.
                For Bogoliubov (spin-wave) vectors use
                :math:`\\eta=\\mathrm{diag}(1,\\ldots,1,-1,\\ldots,-1)`; the sign of
                :math:`\\langle v|\\eta|v\\rangle` is kept.
        """
    def sort(self) -> None:
        ...
    def to_file(self, filename: str, entry: str = 'BZMeshQdc', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/grid", where to write inside the file,
              with a default equal to the object Class name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def BrillouinZone(self) -> BrillouinZone:
        ...
    @property
    def bytes_per_point(self) -> int:
        """
            Return the memory required per interpolation point *result* in bytes
        """
    @property
    def holds_triangulation(self) -> bool:
        """
        Whether the triangulation that refinement works on is in memory (see release_triangulation)
        """
    @property
    def invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def normalizes_vectors(self) -> bool:
        """
            Whether interpolated eigenvectors are normalized, given the stored data; see :py:meth:`set_vector_normalization`
        """
    @property
    def refinable(self) -> bool:
        """
        Whether the mesh can be refined; a mesh read from a file written before refinement existed can't be
        """
    @property
    def refinement_limited(self) -> bool:
        """
        Whether max_points made the mesh coarser than max_size asked for
        """
    @property
    def rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def tetrahedra(self) -> numpy.ndarray[numpy.uint32]:
        ...
    @property
    def values(self) -> numpy.ndarray[numpy.float64]:
        """
            Return a shared view of the stored eigenvalues
        """
    @property
    def vector_metric(self) -> list[float]:
        """
            The diagonal metric used to normalize eigenvectors; empty for the identity
        """
    @property
    def vector_normalization(self) -> str:
        """
            When eigenvectors are normalized: ``"automatic"`` (the default), ``"on"`` or ``"off"``
        """
    @property
    def vectors(self) -> numpy.ndarray[numpy.complex128]:
        """
            Return a shared view of the stored eigenvectors
        """
class BZMeshQcc:
    @staticmethod
    def from_file(filename: str, entry: str = 'BZMeshQcc') -> BZMeshQcc:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/grid", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    def __init__(self, brillouin_zone: BrillouinZone, max_size: float = -1.0, num_levels: int = 3, max_points: int = -1) -> None:
        """
        A structured tetrahedral mesh of a Brillouin zone's irreducible part
        
        A grid of the reciprocal lattice, divided finely enough for ``max_size``, is
        clipped exactly to the irreducible zone.
        
        Parameters
        ----------
        brillouin_zone : BrillouinZone
            The zone whose irreducible part the mesh fills.
        max_size : float, optional (default: -1)
            The largest tetrahedron volume, in cubic reciprocal Angstrom, which sets the
            grid spacing; if not positive, the grid is the reciprocal lattice itself.
            Each grid cell holds six tetrahedra, so ``max_size = node_volume_fraction /
            6`` gives about as many vertices as a :py:class:`BZTrellisQdc` with that
            ``node_volume_fraction``, and ``brillouin_zone.ir_polyhedron.volume / (6 *
            points)`` gives roughly 1.5 to 3 times ``points`` vertices, the most for
            small meshes.
        num_levels : int, optional
            Unused; kept for compatibility.
        max_points : int, optional (default: -1)
            If positive, the grid is coarsened until its estimated number of vertices is
            at most this, with a RuntimeWarning (see :py:attr:`refinement_limited`).
        """
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def __repr__(self) -> str:
        ...
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.complex128], values_elements: numpy.ndarray[numpy.int32], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], sort: bool = False) -> None:
        """
        Provide data required for interpolation to the grid without cost information.
        
        .. Note
        .. ----
        .. This method should probably be followed by :py:meth:`set_cost_info` prior to
        .. any attempt to interpolate the data in the grid.
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below).
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below).
        
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            ``values_data``
        vectors_elements:
            Like ``values_elements`` but for the eigenvectors
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.complex128], values_elements: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
        Provide all data required for interpolation to the grid at once
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below),
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below)
            * which scalar cost function should be used (see below),
            * which vector cost function should be used (see below).
        
            See the note below for the meaning of the last three values.
        values_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvalue elements stored in the grid
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            **values_data**
        vectors_elements:
            Like **values_elements** but for the eigenvectors
        vectors_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvector elements stored in the grid
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        
          Mapping of integers to scalar cost function:
        
          | value | function(x,y) |
          |---|---|
          | 0 | magnitude(x-y) |
        
          Mapping of integers to vector cost function:
        
          | value | function(vec_x, vec_y) |
          |---|---|
          | 0 | sin(hermitian_angle(vec_x, vec_y)) |
          | 1 | vector_distance(vec_x, vec_y) |
          | 2 | 1 - vector_product(vec_x, vec_y) |
          | 3 | vector_angle(vec_x, vec_y) |
          | 4 | hermitian_angle(vec_x, vec_y) |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def ir_interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.complex128], numpy.ndarray[numpy.complex128]]:
        """
          Perform linear interpolation of the stored data at irreducible equivalent points
        
          The irreducible first Brillouin zone is the part of reciprocal space which is
          invariant under application of the integer translations *and* the pointgroup
          operations of a reciprocal space lattice. This method finds points equivalent
          to the input within the irreducible first Brillouin zone and then interpolates
          pre-stored information to provide an estimate at the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              irreducible first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def refine(self, where: typing.Any = None, values: typing.Any = None, vectors: typing.Any = None, resolution: float | None = None, points_per_resolution: float = 2.0) -> numpy.ndarray[numpy.float64]:
        """
        Refine the mesh by bisecting tetrahedra; existing vertices keep their indices and data.
        
        Parameters
        ----------
        where, resolution, points_per_resolution
            As for :py:meth:`refinement_points`, which gives the points this adds.
        values, vectors : numpy.ndarray, optional
            If the mesh holds data (after :py:meth:`fill`), the data for the new vertices,
            laid out per point as the filled data and for exactly the points
            :py:meth:`refinement_points` returns, in that order. Not allowed before the
            mesh is filled.
        
        Returns
        -------
        numpy.ndarray
            The new vertices, shape (N, 3), in relative lattice units, appended to
            :py:attr:`rlu` in this order.
        
        Note
        ----
        The mode permutations found by :py:meth:`sort` are reset; sort again after
        refining if needed.
        """
    def refinement_points(self, where: typing.Any = None, resolution: float | None = None, points_per_resolution: float = 2.0) -> numpy.ndarray[numpy.float64]:
        """
        The points that :py:meth:`refine` would add, without changing the mesh.
        
        Evaluate your model at these points and compare with
        :py:meth:`ir_interpolate_at` there to decide whether refining is worth it; then
        pass the model's values to :py:meth:`refine` with the same arguments.
        
        Parameters
        ----------
        where : None, bool array or int array
            The tetrahedra to split: all of them (None), those where a boolean mask with
            one entry per tetrahedron is true, or those with the given indices.
            Neighbouring tetrahedra are split as needed to keep the mesh conforming, and
            split edges on the zone boundary are split with their symmetry equivalents, so
            that equivalent zone faces keep matching.
        resolution : float, optional
            The resolution limit, in inverse Angstrom. No edge is split to below
            ``resolution / points_per_resolution``: a tetrahedron whose longest edge is at
            most twice that is left whole.
        points_per_resolution : float, optional (default: 2)
            How finely to resolve ``resolution``.
        
        Returns
        -------
        numpy.ndarray
            The new vertices, shape (N, 3), in relative lattice units like :py:attr:`rlu`.
            :py:meth:`refine` appends them to the vertices in this order.
        """
    def release_triangulation(self) -> None:
        """
        Free the memory refinement holds between refinements.
        
        After :py:meth:`refine` (or :py:meth:`refinement_points`) the mesh keeps the
        triangulation refinement works on, several times the memory of the mesh itself.
        This frees it. The mesh is unchanged and can still be refined: the triangulation
        is rebuilt when next needed, which costs a build of the mesh plus a replay of the
        refinements made so far.
        """
    def set_flags_weights(self, values_flags: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_flags: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
          Set :py:class:`~brille._brille.RotatesLike`, :py:class:`~brille._brille.LengthUnit`
          and cost functions plus relative cost weights for the values and vectors
          stored in the object
        
          Parameters
          ----------
          values_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvalues at
              neighbouring grid points for scalar- and vector-like eigenvalues.
          values_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvalue elements stored in the grid
          vectors_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvectors at
              neighbouring grid points for scalar- and vector-like eigenvectors.
          vectors_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvector elements stored in the grid
          sort : bool, optional
              Whether the equivalent-mode permutations should be (re)determined following
              the update to the flags and weights.
        
        
          Note
          ----
            Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:
        
            | value | :py:class:`RotatesLike` |
            |---|---|
            | 0 | `vector` |
            | 1 | `pseudovector` |
            | 2 | `Gamma` |
          
            Mapping of integers to :py:class:`LengthUnit` values:
        
            | value | :py:class:`LengthUnit` |
            |---|---|
            | 0 | `none` |
            | 1 | `angstrom` |
            | 2 | `inverse_angstrom` |
            | 3 | `real_lattice` |
            | 4 | `reciprocal_lattice` |
        
            Mapping of integers to scalar cost function:
          
            | value | function(x,y) |
            |---|---|
            | 0 | magnitude(x-y) |
          
            Mapping of integers to vector cost function:
          
            | value | function(vec_x, vec_y) |
            |---|---|
            | 0 | sin(hermitian_angle(vec_x, vec_y)) |
            | 1 | vector_distance(vec_x, vec_y) |
            | 2 | 1 - vector_product(vec_x, vec_y) |
            | 3 | vector_angle(vec_x, vec_y) |
            | 4 | hermitian_angle(vec_x, vec_y) |
        
            Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def set_vector_normalization(self, normalize: bool | None = True, metric: list[float] | None = None) -> None:
        """
            Choose when interpolated eigenvectors are scaled to unit norm
        
            Linear interpolation between unit eigenvectors gives vectors shorter than
            one wherever neighbouring eigenvectors differ, so structure factors computed
            from them come out too small. Normalization scales each interpolated branch
            :math:`v` to :math:`v/\\sqrt{|\\langle v|M|v\\rangle|}`.
        
            By default it is automatic: eigenvectors stored in Cartesian units
            (:py:class:`LengthUnit` ``angstrom`` or ``inverse_angstrom``, as Euphonic
            stores them) are normalized, and those in lattice units, whose length depends
            on the lattice, are not. The choice survives :py:meth:`fill` and saving to HDF5.
        
            Parameters
            ----------
            normalize : bool or None, optional
                ``True`` always normalizes, and raises a RuntimeError for eigenvectors in
                lattice units; ``False`` never does; ``None`` restores the automatic default.
            metric : float, vector-like, optional
                A diagonal metric :math:`M`, one weight per element of a branch (for
                phonons, :math:`3N`). The default is the identity, the ordinary norm.
                For Bogoliubov (spin-wave) vectors use
                :math:`\\eta=\\mathrm{diag}(1,\\ldots,1,-1,\\ldots,-1)`; the sign of
                :math:`\\langle v|\\eta|v\\rangle` is kept.
        """
    def sort(self) -> None:
        ...
    def to_file(self, filename: str, entry: str = 'BZMeshQcc', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/grid", where to write inside the file,
              with a default equal to the object Class name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def BrillouinZone(self) -> BrillouinZone:
        ...
    @property
    def bytes_per_point(self) -> int:
        """
            Return the memory required per interpolation point *result* in bytes
        """
    @property
    def holds_triangulation(self) -> bool:
        """
        Whether the triangulation that refinement works on is in memory (see release_triangulation)
        """
    @property
    def invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def normalizes_vectors(self) -> bool:
        """
            Whether interpolated eigenvectors are normalized, given the stored data; see :py:meth:`set_vector_normalization`
        """
    @property
    def refinable(self) -> bool:
        """
        Whether the mesh can be refined; a mesh read from a file written before refinement existed can't be
        """
    @property
    def refinement_limited(self) -> bool:
        """
        Whether max_points made the mesh coarser than max_size asked for
        """
    @property
    def rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def tetrahedra(self) -> numpy.ndarray[numpy.uint32]:
        ...
    @property
    def values(self) -> numpy.ndarray[numpy.complex128]:
        """
            Return a shared view of the stored eigenvalues
        """
    @property
    def vector_metric(self) -> list[float]:
        """
            The diagonal metric used to normalize eigenvectors; empty for the identity
        """
    @property
    def vector_normalization(self) -> str:
        """
            When eigenvectors are normalized: ``"automatic"`` (the default), ``"on"`` or ``"off"``
        """
    @property
    def vectors(self) -> numpy.ndarray[numpy.complex128]:
        """
            Return a shared view of the stored eigenvectors
        """
class BZTrellisQdd:
    @staticmethod
    def from_file(filename: str, entry: str = 'BZTrellisQdd') -> BZTrellisQdd:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/grid", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, node_volume_fraction: float = 0.1, always_triangulate: bool = False) -> None:
        """
        A trellis of cubic nodes over a Brillouin zone's irreducible part
        
        Parameters
        ----------
        brillouin_zone : BrillouinZone
            The zone whose irreducible part the trellis fills.
        node_volume_fraction : float, optional (default: 0.1)
            Despite its name, a volume in cubic reciprocal Angstrom, not a fraction:
            the volume of one cubic node, which sets the trellis spacing. To size the
            trellis by its number of points use ``brillouin_zone.ir_polyhedron.volume /
            points``, which gives roughly 1.3 to 2 times ``points`` vertices.
        always_triangulate : bool, optional (default: False)
            Divide every node into tetrahedra, not only those the zone boundary cuts.
        """
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, node_volume_fraction: float, always_triangulate: bool, approx_config: ApproxConfig) -> None:
        ...
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def all_node_types(self) -> list[NodeType]:
        ...
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], vectors_data: numpy.ndarray[numpy.float64], vectors_elements: numpy.ndarray[numpy.int32], sort: bool = False) -> None:
        """
        Provide data required for interpolation to the grid without cost information.
        
        .. Note
        .. ----
        .. This method should probably be followed by :py:meth:`set_cost_info` prior to
        .. any attempt to interpolate the data in the grid.
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below).
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below).
        
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            ``values_data``
        vectors_elements:
            Like ``values_elements`` but for the eigenvectors
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_data: numpy.ndarray[numpy.float64], vectors_elements: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
        Provide all data required for interpolation to the grid at once
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below),
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below)
            * which scalar cost function should be used (see below),
            * which vector cost function should be used (see below).
        
            See the note below for the meaning of the last three values.
        values_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvalue elements stored in the grid
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            **values_data**
        vectors_elements:
            Like **values_elements** but for the eigenvectors
        vectors_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvector elements stored in the grid
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        
          Mapping of integers to scalar cost function:
        
          | value | function(x,y) |
          |---|---|
          | 0 | magnitude(x-y) |
        
          Mapping of integers to vector cost function:
        
          | value | function(vec_x, vec_y) |
          |---|---|
          | 0 | sin(hermitian_angle(vec_x, vec_y)) |
          | 1 | vector_distance(vec_x, vec_y) |
          | 2 | 1 - vector_product(vec_x, vec_y) |
          | 3 | vector_angle(vec_x, vec_y) |
          | 4 | hermitian_angle(vec_x, vec_y) |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.float64], numpy.ndarray[numpy.float64]]:
        """
          Perform linear interpolation of the stored data at equivalent points
        
          The first Brillouin zone is the part of reciprocal space which is invariant
          under application of the integer translations of a reciprocal space lattice.
          This method finds points equivalent to the input within the first Brillouin
          zone and then interpolates pre-stored information to provide an estimate at
          the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method.
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def ir_interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.float64], numpy.ndarray[numpy.float64]]:
        """
          Perform linear interpolation of the stored data at irreducible equivalent points
        
          The irreducible first Brillouin zone is the part of reciprocal space which is
          invariant under application of the integer translations *and* the pointgroup
          operations of a reciprocal space lattice. This method finds points equivalent
          to the input within the irreducible first Brillouin zone and then interpolates
          pre-stored information to provide an estimate at the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              irreducible first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def node_at(self, subscript: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]) -> Polyhedron:
        ...
    def node_at_type(self, subscript: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]) -> NodeType:
        ...
    def node_containing(self, Q: numpy.ndarray[numpy.float64]) -> Polyhedron:
        ...
    def node_containing_type(self, Q: numpy.ndarray[numpy.float64]) -> NodeType:
        ...
    def set_flags_weights(self, values_flags: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_flags: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
          Set :py:class:`~brille._brille.RotatesLike`, :py:class:`~brille._brille.LengthUnit`
          and cost functions plus relative cost weights for the values and vectors
          stored in the object
        
          Parameters
          ----------
          values_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvalues at
              neighbouring grid points for scalar- and vector-like eigenvalues.
          values_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvalue elements stored in the grid
          vectors_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvectors at
              neighbouring grid points for scalar- and vector-like eigenvectors.
          vectors_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvector elements stored in the grid
          sort : bool, optional
              Whether the equivalent-mode permutations should be (re)determined following
              the update to the flags and weights.
        
        
          Note
          ----
            Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:
        
            | value | :py:class:`RotatesLike` |
            |---|---|
            | 0 | `vector` |
            | 1 | `pseudovector` |
            | 2 | `Gamma` |
          
            Mapping of integers to :py:class:`LengthUnit` values:
        
            | value | :py:class:`LengthUnit` |
            |---|---|
            | 0 | `none` |
            | 1 | `angstrom` |
            | 2 | `inverse_angstrom` |
            | 3 | `real_lattice` |
            | 4 | `reciprocal_lattice` |
        
            Mapping of integers to scalar cost function:
          
            | value | function(x,y) |
            |---|---|
            | 0 | magnitude(x-y) |
          
            Mapping of integers to vector cost function:
          
            | value | function(vec_x, vec_y) |
            |---|---|
            | 0 | sin(hermitian_angle(vec_x, vec_y)) |
            | 1 | vector_distance(vec_x, vec_y) |
            | 2 | 1 - vector_product(vec_x, vec_y) |
            | 3 | vector_angle(vec_x, vec_y) |
            | 4 | hermitian_angle(vec_x, vec_y) |
        
            Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def set_vector_normalization(self, normalize: bool | None = True, metric: list[float] | None = None) -> None:
        """
            Choose when interpolated eigenvectors are scaled to unit norm
        
            Linear interpolation between unit eigenvectors gives vectors shorter than
            one wherever neighbouring eigenvectors differ, so structure factors computed
            from them come out too small. Normalization scales each interpolated branch
            :math:`v` to :math:`v/\\sqrt{|\\langle v|M|v\\rangle|}`.
        
            By default it is automatic: eigenvectors stored in Cartesian units
            (:py:class:`LengthUnit` ``angstrom`` or ``inverse_angstrom``, as Euphonic
            stores them) are normalized, and those in lattice units, whose length depends
            on the lattice, are not. The choice survives :py:meth:`fill` and saving to HDF5.
        
            Parameters
            ----------
            normalize : bool or None, optional
                ``True`` always normalizes, and raises a RuntimeError for eigenvectors in
                lattice units; ``False`` never does; ``None`` restores the automatic default.
            metric : float, vector-like, optional
                A diagonal metric :math:`M`, one weight per element of a branch (for
                phonons, :math:`3N`). The default is the identity, the ordinary norm.
                For Bogoliubov (spin-wave) vectors use
                :math:`\\eta=\\mathrm{diag}(1,\\ldots,1,-1,\\ldots,-1)`; the sign of
                :math:`\\langle v|\\eta|v\\rangle` is kept.
        """
    def sort(self) -> None:
        ...
    def to_file(self, filename: str, entry: str = 'BZTrellisQdd', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/grid", where to write inside the file,
              with a default equal to the object Class name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def BrillouinZone(self) -> BrillouinZone:
        ...
    @property
    def bytes_per_point(self) -> int:
        """
            Return the memory required per interpolation point *result* in bytes
        """
    @property
    def inner_invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def inner_rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def normalizes_vectors(self) -> bool:
        """
            Whether interpolated eigenvectors are normalized, given the stored data; see :py:meth:`set_vector_normalization`
        """
    @property
    def outer_invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def outer_rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def tetrahedra(self) -> list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]]:
        ...
    @property
    def values(self) -> numpy.ndarray[numpy.float64]:
        """
            Return a shared view of the stored eigenvalues
        """
    @property
    def vector_metric(self) -> list[float]:
        """
            The diagonal metric used to normalize eigenvectors; empty for the identity
        """
    @property
    def vector_normalization(self) -> str:
        """
            When eigenvectors are normalized: ``"automatic"`` (the default), ``"on"`` or ``"off"``
        """
    @property
    def vectors(self) -> numpy.ndarray[numpy.float64]:
        """
            Return a shared view of the stored eigenvectors
        """
class BZTrellisQdc:
    @staticmethod
    def from_file(filename: str, entry: str = 'BZTrellisQdc') -> BZTrellisQdc:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/grid", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, node_volume_fraction: float = 0.1, always_triangulate: bool = False) -> None:
        """
        A trellis of cubic nodes over a Brillouin zone's irreducible part
        
        Parameters
        ----------
        brillouin_zone : BrillouinZone
            The zone whose irreducible part the trellis fills.
        node_volume_fraction : float, optional (default: 0.1)
            Despite its name, a volume in cubic reciprocal Angstrom, not a fraction:
            the volume of one cubic node, which sets the trellis spacing. To size the
            trellis by its number of points use ``brillouin_zone.ir_polyhedron.volume /
            points``, which gives roughly 1.3 to 2 times ``points`` vertices.
        always_triangulate : bool, optional (default: False)
            Divide every node into tetrahedra, not only those the zone boundary cuts.
        """
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, node_volume_fraction: float, always_triangulate: bool, approx_config: ApproxConfig) -> None:
        ...
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def all_node_types(self) -> list[NodeType]:
        ...
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], sort: bool = False) -> None:
        """
        Provide data required for interpolation to the grid without cost information.
        
        .. Note
        .. ----
        .. This method should probably be followed by :py:meth:`set_cost_info` prior to
        .. any attempt to interpolate the data in the grid.
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below).
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below).
        
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            ``values_data``
        vectors_elements:
            Like ``values_elements`` but for the eigenvectors
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
        Provide all data required for interpolation to the grid at once
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below),
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below)
            * which scalar cost function should be used (see below),
            * which vector cost function should be used (see below).
        
            See the note below for the meaning of the last three values.
        values_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvalue elements stored in the grid
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            **values_data**
        vectors_elements:
            Like **values_elements** but for the eigenvectors
        vectors_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvector elements stored in the grid
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        
          Mapping of integers to scalar cost function:
        
          | value | function(x,y) |
          |---|---|
          | 0 | magnitude(x-y) |
        
          Mapping of integers to vector cost function:
        
          | value | function(vec_x, vec_y) |
          |---|---|
          | 0 | sin(hermitian_angle(vec_x, vec_y)) |
          | 1 | vector_distance(vec_x, vec_y) |
          | 2 | 1 - vector_product(vec_x, vec_y) |
          | 3 | vector_angle(vec_x, vec_y) |
          | 4 | hermitian_angle(vec_x, vec_y) |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.float64], numpy.ndarray[numpy.complex128]]:
        """
          Perform linear interpolation of the stored data at equivalent points
        
          The first Brillouin zone is the part of reciprocal space which is invariant
          under application of the integer translations of a reciprocal space lattice.
          This method finds points equivalent to the input within the first Brillouin
          zone and then interpolates pre-stored information to provide an estimate at
          the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method.
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def ir_interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.float64], numpy.ndarray[numpy.complex128]]:
        """
          Perform linear interpolation of the stored data at irreducible equivalent points
        
          The irreducible first Brillouin zone is the part of reciprocal space which is
          invariant under application of the integer translations *and* the pointgroup
          operations of a reciprocal space lattice. This method finds points equivalent
          to the input within the irreducible first Brillouin zone and then interpolates
          pre-stored information to provide an estimate at the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              irreducible first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def node_at(self, subscript: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]) -> Polyhedron:
        ...
    def node_at_type(self, subscript: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]) -> NodeType:
        ...
    def node_containing(self, Q: numpy.ndarray[numpy.float64]) -> Polyhedron:
        ...
    def node_containing_type(self, Q: numpy.ndarray[numpy.float64]) -> NodeType:
        ...
    def set_flags_weights(self, values_flags: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_flags: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
          Set :py:class:`~brille._brille.RotatesLike`, :py:class:`~brille._brille.LengthUnit`
          and cost functions plus relative cost weights for the values and vectors
          stored in the object
        
          Parameters
          ----------
          values_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvalues at
              neighbouring grid points for scalar- and vector-like eigenvalues.
          values_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvalue elements stored in the grid
          vectors_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvectors at
              neighbouring grid points for scalar- and vector-like eigenvectors.
          vectors_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvector elements stored in the grid
          sort : bool, optional
              Whether the equivalent-mode permutations should be (re)determined following
              the update to the flags and weights.
        
        
          Note
          ----
            Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:
        
            | value | :py:class:`RotatesLike` |
            |---|---|
            | 0 | `vector` |
            | 1 | `pseudovector` |
            | 2 | `Gamma` |
          
            Mapping of integers to :py:class:`LengthUnit` values:
        
            | value | :py:class:`LengthUnit` |
            |---|---|
            | 0 | `none` |
            | 1 | `angstrom` |
            | 2 | `inverse_angstrom` |
            | 3 | `real_lattice` |
            | 4 | `reciprocal_lattice` |
        
            Mapping of integers to scalar cost function:
          
            | value | function(x,y) |
            |---|---|
            | 0 | magnitude(x-y) |
          
            Mapping of integers to vector cost function:
          
            | value | function(vec_x, vec_y) |
            |---|---|
            | 0 | sin(hermitian_angle(vec_x, vec_y)) |
            | 1 | vector_distance(vec_x, vec_y) |
            | 2 | 1 - vector_product(vec_x, vec_y) |
            | 3 | vector_angle(vec_x, vec_y) |
            | 4 | hermitian_angle(vec_x, vec_y) |
        
            Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def set_vector_normalization(self, normalize: bool | None = True, metric: list[float] | None = None) -> None:
        """
            Choose when interpolated eigenvectors are scaled to unit norm
        
            Linear interpolation between unit eigenvectors gives vectors shorter than
            one wherever neighbouring eigenvectors differ, so structure factors computed
            from them come out too small. Normalization scales each interpolated branch
            :math:`v` to :math:`v/\\sqrt{|\\langle v|M|v\\rangle|}`.
        
            By default it is automatic: eigenvectors stored in Cartesian units
            (:py:class:`LengthUnit` ``angstrom`` or ``inverse_angstrom``, as Euphonic
            stores them) are normalized, and those in lattice units, whose length depends
            on the lattice, are not. The choice survives :py:meth:`fill` and saving to HDF5.
        
            Parameters
            ----------
            normalize : bool or None, optional
                ``True`` always normalizes, and raises a RuntimeError for eigenvectors in
                lattice units; ``False`` never does; ``None`` restores the automatic default.
            metric : float, vector-like, optional
                A diagonal metric :math:`M`, one weight per element of a branch (for
                phonons, :math:`3N`). The default is the identity, the ordinary norm.
                For Bogoliubov (spin-wave) vectors use
                :math:`\\eta=\\mathrm{diag}(1,\\ldots,1,-1,\\ldots,-1)`; the sign of
                :math:`\\langle v|\\eta|v\\rangle` is kept.
        """
    def sort(self) -> None:
        ...
    def to_file(self, filename: str, entry: str = 'BZTrellisQdc', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/grid", where to write inside the file,
              with a default equal to the object Class name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def BrillouinZone(self) -> BrillouinZone:
        ...
    @property
    def bytes_per_point(self) -> int:
        """
            Return the memory required per interpolation point *result* in bytes
        """
    @property
    def inner_invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def inner_rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def normalizes_vectors(self) -> bool:
        """
            Whether interpolated eigenvectors are normalized, given the stored data; see :py:meth:`set_vector_normalization`
        """
    @property
    def outer_invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def outer_rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def tetrahedra(self) -> list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]]:
        ...
    @property
    def values(self) -> numpy.ndarray[numpy.float64]:
        """
            Return a shared view of the stored eigenvalues
        """
    @property
    def vector_metric(self) -> list[float]:
        """
            The diagonal metric used to normalize eigenvectors; empty for the identity
        """
    @property
    def vector_normalization(self) -> str:
        """
            When eigenvectors are normalized: ``"automatic"`` (the default), ``"on"`` or ``"off"``
        """
    @property
    def vectors(self) -> numpy.ndarray[numpy.complex128]:
        """
            Return a shared view of the stored eigenvectors
        """
class BZTrellisQcc:
    @staticmethod
    def from_file(filename: str, entry: str = 'BZTrellisQcc') -> BZTrellisQcc:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/grid", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, node_volume_fraction: float = 0.1, always_triangulate: bool = False) -> None:
        """
        A trellis of cubic nodes over a Brillouin zone's irreducible part
        
        Parameters
        ----------
        brillouin_zone : BrillouinZone
            The zone whose irreducible part the trellis fills.
        node_volume_fraction : float, optional (default: 0.1)
            Despite its name, a volume in cubic reciprocal Angstrom, not a fraction:
            the volume of one cubic node, which sets the trellis spacing. To size the
            trellis by its number of points use ``brillouin_zone.ir_polyhedron.volume /
            points``, which gives roughly 1.3 to 2 times ``points`` vertices.
        always_triangulate : bool, optional (default: False)
            Divide every node into tetrahedra, not only those the zone boundary cuts.
        """
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, node_volume_fraction: float, always_triangulate: bool, approx_config: ApproxConfig) -> None:
        ...
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def all_node_types(self) -> list[NodeType]:
        ...
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.complex128], values_elements: numpy.ndarray[numpy.int32], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], sort: bool = False) -> None:
        """
        Provide data required for interpolation to the grid without cost information.
        
        .. Note
        .. ----
        .. This method should probably be followed by :py:meth:`set_cost_info` prior to
        .. any attempt to interpolate the data in the grid.
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below).
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below).
        
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            ``values_data``
        vectors_elements:
            Like ``values_elements`` but for the eigenvectors
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.complex128], values_elements: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
        Provide all data required for interpolation to the grid at once
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below),
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below)
            * which scalar cost function should be used (see below),
            * which vector cost function should be used (see below).
        
            See the note below for the meaning of the last three values.
        values_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvalue elements stored in the grid
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            **values_data**
        vectors_elements:
            Like **values_elements** but for the eigenvectors
        vectors_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvector elements stored in the grid
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        
          Mapping of integers to scalar cost function:
        
          | value | function(x,y) |
          |---|---|
          | 0 | magnitude(x-y) |
        
          Mapping of integers to vector cost function:
        
          | value | function(vec_x, vec_y) |
          |---|---|
          | 0 | sin(hermitian_angle(vec_x, vec_y)) |
          | 1 | vector_distance(vec_x, vec_y) |
          | 2 | 1 - vector_product(vec_x, vec_y) |
          | 3 | vector_angle(vec_x, vec_y) |
          | 4 | hermitian_angle(vec_x, vec_y) |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.complex128], numpy.ndarray[numpy.complex128]]:
        """
          Perform linear interpolation of the stored data at equivalent points
        
          The first Brillouin zone is the part of reciprocal space which is invariant
          under application of the integer translations of a reciprocal space lattice.
          This method finds points equivalent to the input within the first Brillouin
          zone and then interpolates pre-stored information to provide an estimate at
          the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method.
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def ir_interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.complex128], numpy.ndarray[numpy.complex128]]:
        """
          Perform linear interpolation of the stored data at irreducible equivalent points
        
          The irreducible first Brillouin zone is the part of reciprocal space which is
          invariant under application of the integer translations *and* the pointgroup
          operations of a reciprocal space lattice. This method finds points equivalent
          to the input within the irreducible first Brillouin zone and then interpolates
          pre-stored information to provide an estimate at the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              irreducible first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def node_at(self, subscript: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]) -> Polyhedron:
        ...
    def node_at_type(self, subscript: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]) -> NodeType:
        ...
    def node_containing(self, Q: numpy.ndarray[numpy.float64]) -> Polyhedron:
        ...
    def node_containing_type(self, Q: numpy.ndarray[numpy.float64]) -> NodeType:
        ...
    def set_flags_weights(self, values_flags: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_flags: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
          Set :py:class:`~brille._brille.RotatesLike`, :py:class:`~brille._brille.LengthUnit`
          and cost functions plus relative cost weights for the values and vectors
          stored in the object
        
          Parameters
          ----------
          values_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvalues at
              neighbouring grid points for scalar- and vector-like eigenvalues.
          values_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvalue elements stored in the grid
          vectors_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvectors at
              neighbouring grid points for scalar- and vector-like eigenvectors.
          vectors_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvector elements stored in the grid
          sort : bool, optional
              Whether the equivalent-mode permutations should be (re)determined following
              the update to the flags and weights.
        
        
          Note
          ----
            Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:
        
            | value | :py:class:`RotatesLike` |
            |---|---|
            | 0 | `vector` |
            | 1 | `pseudovector` |
            | 2 | `Gamma` |
          
            Mapping of integers to :py:class:`LengthUnit` values:
        
            | value | :py:class:`LengthUnit` |
            |---|---|
            | 0 | `none` |
            | 1 | `angstrom` |
            | 2 | `inverse_angstrom` |
            | 3 | `real_lattice` |
            | 4 | `reciprocal_lattice` |
        
            Mapping of integers to scalar cost function:
          
            | value | function(x,y) |
            |---|---|
            | 0 | magnitude(x-y) |
          
            Mapping of integers to vector cost function:
          
            | value | function(vec_x, vec_y) |
            |---|---|
            | 0 | sin(hermitian_angle(vec_x, vec_y)) |
            | 1 | vector_distance(vec_x, vec_y) |
            | 2 | 1 - vector_product(vec_x, vec_y) |
            | 3 | vector_angle(vec_x, vec_y) |
            | 4 | hermitian_angle(vec_x, vec_y) |
        
            Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def set_vector_normalization(self, normalize: bool | None = True, metric: list[float] | None = None) -> None:
        """
            Choose when interpolated eigenvectors are scaled to unit norm
        
            Linear interpolation between unit eigenvectors gives vectors shorter than
            one wherever neighbouring eigenvectors differ, so structure factors computed
            from them come out too small. Normalization scales each interpolated branch
            :math:`v` to :math:`v/\\sqrt{|\\langle v|M|v\\rangle|}`.
        
            By default it is automatic: eigenvectors stored in Cartesian units
            (:py:class:`LengthUnit` ``angstrom`` or ``inverse_angstrom``, as Euphonic
            stores them) are normalized, and those in lattice units, whose length depends
            on the lattice, are not. The choice survives :py:meth:`fill` and saving to HDF5.
        
            Parameters
            ----------
            normalize : bool or None, optional
                ``True`` always normalizes, and raises a RuntimeError for eigenvectors in
                lattice units; ``False`` never does; ``None`` restores the automatic default.
            metric : float, vector-like, optional
                A diagonal metric :math:`M`, one weight per element of a branch (for
                phonons, :math:`3N`). The default is the identity, the ordinary norm.
                For Bogoliubov (spin-wave) vectors use
                :math:`\\eta=\\mathrm{diag}(1,\\ldots,1,-1,\\ldots,-1)`; the sign of
                :math:`\\langle v|\\eta|v\\rangle` is kept.
        """
    def sort(self) -> None:
        ...
    def to_file(self, filename: str, entry: str = 'BZTrellisQcc', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/grid", where to write inside the file,
              with a default equal to the object Class name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def BrillouinZone(self) -> BrillouinZone:
        ...
    @property
    def bytes_per_point(self) -> int:
        """
            Return the memory required per interpolation point *result* in bytes
        """
    @property
    def inner_invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def inner_rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def normalizes_vectors(self) -> bool:
        """
            Whether interpolated eigenvectors are normalized, given the stored data; see :py:meth:`set_vector_normalization`
        """
    @property
    def outer_invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def outer_rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def tetrahedra(self) -> list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]]:
        ...
    @property
    def values(self) -> numpy.ndarray[numpy.complex128]:
        """
            Return a shared view of the stored eigenvalues
        """
    @property
    def vector_metric(self) -> list[float]:
        """
            The diagonal metric used to normalize eigenvectors; empty for the identity
        """
    @property
    def vector_normalization(self) -> str:
        """
            When eigenvectors are normalized: ``"automatic"`` (the default), ``"on"`` or ``"off"``
        """
    @property
    def vectors(self) -> numpy.ndarray[numpy.complex128]:
        """
            Return a shared view of the stored eigenvectors
        """
class BZNestQdd:
    @staticmethod
    def from_file(filename: str, entry: str = 'BZNestQdd') -> BZNestQdd:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/grid", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, max_volume: float, max_branchings: int = 5) -> None:
        ...
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, number_density: int, max_branchings: int = 5) -> None:
        ...
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], vectors_data: numpy.ndarray[numpy.float64], vectors_elements: numpy.ndarray[numpy.int32], sort: bool = False) -> None:
        """
        Provide data required for interpolation to the grid without cost information.
        
        .. Note
        .. ----
        .. This method should probably be followed by :py:meth:`set_cost_info` prior to
        .. any attempt to interpolate the data in the grid.
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below).
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below).
        
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            ``values_data``
        vectors_elements:
            Like ``values_elements`` but for the eigenvectors
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_data: numpy.ndarray[numpy.float64], vectors_elements: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
        Provide all data required for interpolation to the grid at once
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below),
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below)
            * which scalar cost function should be used (see below),
            * which vector cost function should be used (see below).
        
            See the note below for the meaning of the last three values.
        values_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvalue elements stored in the grid
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            **values_data**
        vectors_elements:
            Like **values_elements** but for the eigenvectors
        vectors_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvector elements stored in the grid
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        
          Mapping of integers to scalar cost function:
        
          | value | function(x,y) |
          |---|---|
          | 0 | magnitude(x-y) |
        
          Mapping of integers to vector cost function:
        
          | value | function(vec_x, vec_y) |
          |---|---|
          | 0 | sin(hermitian_angle(vec_x, vec_y)) |
          | 1 | vector_distance(vec_x, vec_y) |
          | 2 | 1 - vector_product(vec_x, vec_y) |
          | 3 | vector_angle(vec_x, vec_y) |
          | 4 | hermitian_angle(vec_x, vec_y) |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def ir_interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.float64], numpy.ndarray[numpy.float64]]:
        """
          Perform linear interpolation of the stored data at irreducible equivalent points
        
          The irreducible first Brillouin zone is the part of reciprocal space which is
          invariant under application of the integer translations *and* the pointgroup
          operations of a reciprocal space lattice. This method finds points equivalent
          to the input within the irreducible first Brillouin zone and then interpolates
          pre-stored information to provide an estimate at the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              irreducible first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def set_flags_weights(self, values_flags: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_flags: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
          Set :py:class:`~brille._brille.RotatesLike`, :py:class:`~brille._brille.LengthUnit`
          and cost functions plus relative cost weights for the values and vectors
          stored in the object
        
          Parameters
          ----------
          values_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvalues at
              neighbouring grid points for scalar- and vector-like eigenvalues.
          values_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvalue elements stored in the grid
          vectors_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvectors at
              neighbouring grid points for scalar- and vector-like eigenvectors.
          vectors_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvector elements stored in the grid
          sort : bool, optional
              Whether the equivalent-mode permutations should be (re)determined following
              the update to the flags and weights.
        
        
          Note
          ----
            Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:
        
            | value | :py:class:`RotatesLike` |
            |---|---|
            | 0 | `vector` |
            | 1 | `pseudovector` |
            | 2 | `Gamma` |
          
            Mapping of integers to :py:class:`LengthUnit` values:
        
            | value | :py:class:`LengthUnit` |
            |---|---|
            | 0 | `none` |
            | 1 | `angstrom` |
            | 2 | `inverse_angstrom` |
            | 3 | `real_lattice` |
            | 4 | `reciprocal_lattice` |
        
            Mapping of integers to scalar cost function:
          
            | value | function(x,y) |
            |---|---|
            | 0 | magnitude(x-y) |
          
            Mapping of integers to vector cost function:
          
            | value | function(vec_x, vec_y) |
            |---|---|
            | 0 | sin(hermitian_angle(vec_x, vec_y)) |
            | 1 | vector_distance(vec_x, vec_y) |
            | 2 | 1 - vector_product(vec_x, vec_y) |
            | 3 | vector_angle(vec_x, vec_y) |
            | 4 | hermitian_angle(vec_x, vec_y) |
        
            Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def set_vector_normalization(self, normalize: bool | None = True, metric: list[float] | None = None) -> None:
        """
            Choose when interpolated eigenvectors are scaled to unit norm
        
            Linear interpolation between unit eigenvectors gives vectors shorter than
            one wherever neighbouring eigenvectors differ, so structure factors computed
            from them come out too small. Normalization scales each interpolated branch
            :math:`v` to :math:`v/\\sqrt{|\\langle v|M|v\\rangle|}`.
        
            By default it is automatic: eigenvectors stored in Cartesian units
            (:py:class:`LengthUnit` ``angstrom`` or ``inverse_angstrom``, as Euphonic
            stores them) are normalized, and those in lattice units, whose length depends
            on the lattice, are not. The choice survives :py:meth:`fill` and saving to HDF5.
        
            Parameters
            ----------
            normalize : bool or None, optional
                ``True`` always normalizes, and raises a RuntimeError for eigenvectors in
                lattice units; ``False`` never does; ``None`` restores the automatic default.
            metric : float, vector-like, optional
                A diagonal metric :math:`M`, one weight per element of a branch (for
                phonons, :math:`3N`). The default is the identity, the ordinary norm.
                For Bogoliubov (spin-wave) vectors use
                :math:`\\eta=\\mathrm{diag}(1,\\ldots,1,-1,\\ldots,-1)`; the sign of
                :math:`\\langle v|\\eta|v\\rangle` is kept.
        """
    def sort(self) -> None:
        ...
    def to_file(self, filename: str, entry: str = 'BZNestQdd', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/grid", where to write inside the file,
              with a default equal to the object Class name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def BrillouinZone(self) -> BrillouinZone:
        ...
    @property
    def all_invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def all_rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def bytes_per_point(self) -> int:
        """
            Return the memory required per interpolation point *result* in bytes
        """
    @property
    def invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def normalizes_vectors(self) -> bool:
        """
            Whether interpolated eigenvectors are normalized, given the stored data; see :py:meth:`set_vector_normalization`
        """
    @property
    def rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def tetrahedra(self) -> list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]]:
        ...
    @property
    def values(self) -> numpy.ndarray[numpy.float64]:
        """
            Return a shared view of the stored eigenvalues
        """
    @property
    def vector_metric(self) -> list[float]:
        """
            The diagonal metric used to normalize eigenvectors; empty for the identity
        """
    @property
    def vector_normalization(self) -> str:
        """
            When eigenvectors are normalized: ``"automatic"`` (the default), ``"on"`` or ``"off"``
        """
    @property
    def vectors(self) -> numpy.ndarray[numpy.float64]:
        """
            Return a shared view of the stored eigenvectors
        """
class BZNestQdc:
    @staticmethod
    def from_file(filename: str, entry: str = 'BZNestQdc') -> BZNestQdc:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/grid", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, max_volume: float, max_branchings: int = 5) -> None:
        ...
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, number_density: int, max_branchings: int = 5) -> None:
        ...
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], sort: bool = False) -> None:
        """
        Provide data required for interpolation to the grid without cost information.
        
        .. Note
        .. ----
        .. This method should probably be followed by :py:meth:`set_cost_info` prior to
        .. any attempt to interpolate the data in the grid.
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below).
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below).
        
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            ``values_data``
        vectors_elements:
            Like ``values_elements`` but for the eigenvectors
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.float64], values_elements: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
        Provide all data required for interpolation to the grid at once
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below),
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below)
            * which scalar cost function should be used (see below),
            * which vector cost function should be used (see below).
        
            See the note below for the meaning of the last three values.
        values_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvalue elements stored in the grid
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            **values_data**
        vectors_elements:
            Like **values_elements** but for the eigenvectors
        vectors_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvector elements stored in the grid
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        
          Mapping of integers to scalar cost function:
        
          | value | function(x,y) |
          |---|---|
          | 0 | magnitude(x-y) |
        
          Mapping of integers to vector cost function:
        
          | value | function(vec_x, vec_y) |
          |---|---|
          | 0 | sin(hermitian_angle(vec_x, vec_y)) |
          | 1 | vector_distance(vec_x, vec_y) |
          | 2 | 1 - vector_product(vec_x, vec_y) |
          | 3 | vector_angle(vec_x, vec_y) |
          | 4 | hermitian_angle(vec_x, vec_y) |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def ir_interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.float64], numpy.ndarray[numpy.complex128]]:
        """
          Perform linear interpolation of the stored data at irreducible equivalent points
        
          The irreducible first Brillouin zone is the part of reciprocal space which is
          invariant under application of the integer translations *and* the pointgroup
          operations of a reciprocal space lattice. This method finds points equivalent
          to the input within the irreducible first Brillouin zone and then interpolates
          pre-stored information to provide an estimate at the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              irreducible first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def set_flags_weights(self, values_flags: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_flags: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
          Set :py:class:`~brille._brille.RotatesLike`, :py:class:`~brille._brille.LengthUnit`
          and cost functions plus relative cost weights for the values and vectors
          stored in the object
        
          Parameters
          ----------
          values_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvalues at
              neighbouring grid points for scalar- and vector-like eigenvalues.
          values_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvalue elements stored in the grid
          vectors_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvectors at
              neighbouring grid points for scalar- and vector-like eigenvectors.
          vectors_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvector elements stored in the grid
          sort : bool, optional
              Whether the equivalent-mode permutations should be (re)determined following
              the update to the flags and weights.
        
        
          Note
          ----
            Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:
        
            | value | :py:class:`RotatesLike` |
            |---|---|
            | 0 | `vector` |
            | 1 | `pseudovector` |
            | 2 | `Gamma` |
          
            Mapping of integers to :py:class:`LengthUnit` values:
        
            | value | :py:class:`LengthUnit` |
            |---|---|
            | 0 | `none` |
            | 1 | `angstrom` |
            | 2 | `inverse_angstrom` |
            | 3 | `real_lattice` |
            | 4 | `reciprocal_lattice` |
        
            Mapping of integers to scalar cost function:
          
            | value | function(x,y) |
            |---|---|
            | 0 | magnitude(x-y) |
          
            Mapping of integers to vector cost function:
          
            | value | function(vec_x, vec_y) |
            |---|---|
            | 0 | sin(hermitian_angle(vec_x, vec_y)) |
            | 1 | vector_distance(vec_x, vec_y) |
            | 2 | 1 - vector_product(vec_x, vec_y) |
            | 3 | vector_angle(vec_x, vec_y) |
            | 4 | hermitian_angle(vec_x, vec_y) |
        
            Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def set_vector_normalization(self, normalize: bool | None = True, metric: list[float] | None = None) -> None:
        """
            Choose when interpolated eigenvectors are scaled to unit norm
        
            Linear interpolation between unit eigenvectors gives vectors shorter than
            one wherever neighbouring eigenvectors differ, so structure factors computed
            from them come out too small. Normalization scales each interpolated branch
            :math:`v` to :math:`v/\\sqrt{|\\langle v|M|v\\rangle|}`.
        
            By default it is automatic: eigenvectors stored in Cartesian units
            (:py:class:`LengthUnit` ``angstrom`` or ``inverse_angstrom``, as Euphonic
            stores them) are normalized, and those in lattice units, whose length depends
            on the lattice, are not. The choice survives :py:meth:`fill` and saving to HDF5.
        
            Parameters
            ----------
            normalize : bool or None, optional
                ``True`` always normalizes, and raises a RuntimeError for eigenvectors in
                lattice units; ``False`` never does; ``None`` restores the automatic default.
            metric : float, vector-like, optional
                A diagonal metric :math:`M`, one weight per element of a branch (for
                phonons, :math:`3N`). The default is the identity, the ordinary norm.
                For Bogoliubov (spin-wave) vectors use
                :math:`\\eta=\\mathrm{diag}(1,\\ldots,1,-1,\\ldots,-1)`; the sign of
                :math:`\\langle v|\\eta|v\\rangle` is kept.
        """
    def sort(self) -> None:
        ...
    def to_file(self, filename: str, entry: str = 'BZNestQdc', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/grid", where to write inside the file,
              with a default equal to the object Class name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def BrillouinZone(self) -> BrillouinZone:
        ...
    @property
    def all_invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def all_rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def bytes_per_point(self) -> int:
        """
            Return the memory required per interpolation point *result* in bytes
        """
    @property
    def invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def normalizes_vectors(self) -> bool:
        """
            Whether interpolated eigenvectors are normalized, given the stored data; see :py:meth:`set_vector_normalization`
        """
    @property
    def rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def tetrahedra(self) -> list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]]:
        ...
    @property
    def values(self) -> numpy.ndarray[numpy.float64]:
        """
            Return a shared view of the stored eigenvalues
        """
    @property
    def vector_metric(self) -> list[float]:
        """
            The diagonal metric used to normalize eigenvectors; empty for the identity
        """
    @property
    def vector_normalization(self) -> str:
        """
            When eigenvectors are normalized: ``"automatic"`` (the default), ``"on"`` or ``"off"``
        """
    @property
    def vectors(self) -> numpy.ndarray[numpy.complex128]:
        """
            Return a shared view of the stored eigenvectors
        """
class BZNestQcc:
    @staticmethod
    def from_file(filename: str, entry: str = 'BZNestQcc') -> BZNestQcc:
        """
          Load an object from an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to read from
          entry: str
              The group path, e.g., "my/cool/grid", where to read from inside the file,
              with a default equal to the object Class name
        
          Returns
          -------
          clsObj
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, max_volume: float, max_branchings: int = 5) -> None:
        ...
    @typing.overload
    def __init__(self, brillouin_zone: BrillouinZone, number_density: int, max_branchings: int = 5) -> None:
        ...
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.complex128], values_elements: numpy.ndarray[numpy.int32], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], sort: bool = False) -> None:
        """
        Provide data required for interpolation to the grid without cost information.
        
        .. Note
        .. ----
        .. This method should probably be followed by :py:meth:`set_cost_info` prior to
        .. any attempt to interpolate the data in the grid.
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below).
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below).
        
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            ``values_data``
        vectors_elements:
            Like ``values_elements`` but for the eigenvectors
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        """
    @typing.overload
    def fill(self, values_data: numpy.ndarray[numpy.complex128], values_elements: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_data: numpy.ndarray[numpy.complex128], vectors_elements: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
        Provide all data required for interpolation to the grid at once
        
        Parameters
        ----------
        values_data : :py:class:`numpy.ndarray`
            The eigenvalue data to be stored in the grid. The first dimension must be
            equal in size to the number of grid-vertices. If two dimensional the second
            dimension is interpreted as all information for a single mode flattened
            and concatenated into (scalars, vectors, matrices) -- in that order.
            If more than two dimensional, the second dimension indexes modes and
            higher dimensions will be flattened *as if row ordered* and must flatten into
            a concatenated list of (scalars, vectors, matrices).
            If the provided array can be interpreted as a contiguous row-ordered two
            dimensional array it will be used in place, otherwise a copy will be made.
        values_elements: integer vector-like
            A multi-purpose vector containing, in order:
        
            * the number of scalar-like eigenvalue elements,
            * the number of vector-like eigenvalue *elements* (must be :math:`3\\times N`),
            * the number of matrix-like eigenvalue *elements* (must be :math:`9\\times N`),
            * an integer :py:class:`RotatesLike` value denoting
            *how* the vector-like and matrix-like parts transform under application
            of a symmetry operation (see note below),
            * an integer :py:class:`LengthUnit` value denoting what units
            the vector-like and matrix-like parts are in (see note below)
            * which scalar cost function should be used (see below),
            * which vector cost function should be used (see below).
        
            See the note below for the meaning of the last three values.
        values_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvalue elements stored in the grid
        vectors_data : :py:class:`numpy.ndarray`
            The eigenvector data to be stored in the grid. Same shape restrictions as
            **values_data**
        vectors_elements:
            Like **values_elements** but for the eigenvectors
        vectors_weights : float, vector-like
            The relative cost weights between scalar-, vector-, and matrix- like
            eigenvector elements stored in the grid
        sort : logical (default ``False``)
            Whether the equivalent-mode permutations should be (re)determined following
            the update to the flags and weights.
        
        
        Note
        ----
          Mapping of integers to :py:class:`RotatesLike` values:
        
          | value | :py:class:`RotatesLike` |
          |---|---|
          | 0 | `vector` |
          | 1 | `pseudovector` |
          | 2 | `Gamma` |
        
          Mapping of integers to :py:class:`LengthUnit` values:
        
          | value | :py:class:`LengthUnit` |
          |---|---|
          | 0 | `none` |
          | 1 | `angstrom` |
          | 2 | `inverse_angstrom` |
          | 3 | `real_lattice` |
          | 4 | `reciprocal_lattice` |
        
          Integer values outside of the mapped range (or missing) are replaced by 3.
        
          Phonon eigenvectors (:py:class:`RotatesLike` `Gamma`) must use the "cell"
          phase convention, in which they are periodic in reciprocal space.
          Eigenvectors in the "atom" convention, which phonopy uses, give wrong results
          without an error; see :ref:`phase_convention` for how to convert them.
        
          Mapping of integers to scalar cost function:
        
          | value | function(x,y) |
          |---|---|
          | 0 | magnitude(x-y) |
        
          Mapping of integers to vector cost function:
        
          | value | function(vec_x, vec_y) |
          |---|---|
          | 0 | sin(hermitian_angle(vec_x, vec_y)) |
          | 1 | vector_distance(vec_x, vec_y) |
          | 2 | 1 - vector_product(vec_x, vec_y) |
          | 3 | vector_angle(vec_x, vec_y) |
          | 4 | hermitian_angle(vec_x, vec_y) |
        
          Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def ir_interpolate_at(self, Q: numpy.ndarray[numpy.float64], useparallel: bool = False, threads: int = -1, do_not_move_points: bool = False) -> tuple[numpy.ndarray[numpy.complex128], numpy.ndarray[numpy.complex128]]:
        """
          Perform linear interpolation of the stored data at irreducible equivalent points
        
          The irreducible first Brillouin zone is the part of reciprocal space which is
          invariant under application of the integer translations *and* the pointgroup
          operations of a reciprocal space lattice. This method finds points equivalent
          to the input within the irreducible first Brillouin zone and then interpolates
          pre-stored information to provide an estimate at the found positions.
        
          Parameters
          ----------
          Q : :py:class:`numpy.ndarray`
              A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
              which an interpolated result is required, expressed in units of the
              reciprocal lattice.
          useparallel : bool, optional
              Whether a serial or parallel code should be utilised
          threads : int, optional
              How many parallel threads should be utilised; if this value is less than one,
              the ``BRILLE_NUM_THREADS`` environment variable sets the number, or one
              thread per logical core is used if it is not set.
          do_not_move_points: bool, optional
              If ``True`` the provided **Q** points must already lie within the first Brillouin
              zone. No check is made to verify this requirement and if any **Q** lie outside
              of the gridded volume out-of-bounds errors may result in bad data or runtime
              errors.
        
          Returns
          -------
          tuple
              The interpolated eigenvalues and eigenvectors at the equivalent
              irreducible first Brillouin zone points.
              The shape of each output will depend on the shape of the data provided to
              the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
              If the filled eigenvalues were of shape
              ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
              ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
              ``[N_Q_points, 3]`` then the output shapes will be
              ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
              for the eigenvalues and eigenvectors, respectively.
        """
    def set_flags_weights(self, values_flags: numpy.ndarray[numpy.int32], values_weights: numpy.ndarray[numpy.float64], vectors_flags: numpy.ndarray[numpy.int32], vectors_weights: numpy.ndarray[numpy.float64], sort: bool = False) -> None:
        """
          Set :py:class:`~brille._brille.RotatesLike`, :py:class:`~brille._brille.LengthUnit`
          and cost functions plus relative cost weights for the values and vectors
          stored in the object
        
          Parameters
          ----------
          values_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvalues at
              neighbouring grid points for scalar- and vector-like eigenvalues.
          values_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvalue elements stored in the grid
          vectors_flags : integer, vector-like
              One or more values indicating the :py:class:`~brille._brille.RotatesLike`
              value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
              value, plus which cost function to use when comparing stored eigenvectors at
              neighbouring grid points for scalar- and vector-like eigenvectors.
          vectors_weights : float, vector-like
              The relative cost weights between scalar-, vector-, and matrix- like
              eigenvector elements stored in the grid
          sort : bool, optional
              Whether the equivalent-mode permutations should be (re)determined following
              the update to the flags and weights.
        
        
          Note
          ----
            Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:
        
            | value | :py:class:`RotatesLike` |
            |---|---|
            | 0 | `vector` |
            | 1 | `pseudovector` |
            | 2 | `Gamma` |
          
            Mapping of integers to :py:class:`LengthUnit` values:
        
            | value | :py:class:`LengthUnit` |
            |---|---|
            | 0 | `none` |
            | 1 | `angstrom` |
            | 2 | `inverse_angstrom` |
            | 3 | `real_lattice` |
            | 4 | `reciprocal_lattice` |
        
            Mapping of integers to scalar cost function:
          
            | value | function(x,y) |
            |---|---|
            | 0 | magnitude(x-y) |
          
            Mapping of integers to vector cost function:
          
            | value | function(vec_x, vec_y) |
            |---|---|
            | 0 | sin(hermitian_angle(vec_x, vec_y)) |
            | 1 | vector_distance(vec_x, vec_y) |
            | 2 | 1 - vector_product(vec_x, vec_y) |
            | 3 | vector_angle(vec_x, vec_y) |
            | 4 | hermitian_angle(vec_x, vec_y) |
        
            Integer values outside of the mapped range (or missing) are replaced by 0.
        """
    def set_vector_normalization(self, normalize: bool | None = True, metric: list[float] | None = None) -> None:
        """
            Choose when interpolated eigenvectors are scaled to unit norm
        
            Linear interpolation between unit eigenvectors gives vectors shorter than
            one wherever neighbouring eigenvectors differ, so structure factors computed
            from them come out too small. Normalization scales each interpolated branch
            :math:`v` to :math:`v/\\sqrt{|\\langle v|M|v\\rangle|}`.
        
            By default it is automatic: eigenvectors stored in Cartesian units
            (:py:class:`LengthUnit` ``angstrom`` or ``inverse_angstrom``, as Euphonic
            stores them) are normalized, and those in lattice units, whose length depends
            on the lattice, are not. The choice survives :py:meth:`fill` and saving to HDF5.
        
            Parameters
            ----------
            normalize : bool or None, optional
                ``True`` always normalizes, and raises a RuntimeError for eigenvectors in
                lattice units; ``False`` never does; ``None`` restores the automatic default.
            metric : float, vector-like, optional
                A diagonal metric :math:`M`, one weight per element of a branch (for
                phonons, :math:`3N`). The default is the identity, the ordinary norm.
                For Bogoliubov (spin-wave) vectors use
                :math:`\\eta=\\mathrm{diag}(1,\\ldots,1,-1,\\ldots,-1)`; the sign of
                :math:`\\langle v|\\eta|v\\rangle` is kept.
        """
    def sort(self) -> None:
        ...
    def to_file(self, filename: str, entry: str = 'BZNestQcc', flags: str = 'ac') -> bool:
        """
          Save the object to an HDF5 file
        
          Parameters
          ----------
          filename : str
              The full path specification for the file to write into
          entry: str
              The group path, e.g., "my/cool/grid", where to write inside the file,
              with a default equal to the object Class name
          flags: str
              The HDF5 permissions to use when opening the file. Default 'a' writes to an
              existing file -- if `entry` exists in the file it is overwritten.
        
          Note
          ----
          Possible `flags` are:
        
          | `flags` | meaning | HDF equivalent |
          |---|---|---|
          | 'r' | read | H5F_ACC_RDONLY |
          | 'x' | write, error if exists | H5F_ACC_EXCL |
          | 'a' | write, append to file | H5F_ACC_RDWR |
          | 'c' | write, error if exists | H5F_ACC_CREAT |
          | 't' | write, replace existing | H5F_ACC_TRUNC |
        
        
          Returns
          -------
          bool
              Indication of writing success.
        """
    @property
    def BrillouinZone(self) -> BrillouinZone:
        ...
    @property
    def all_invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def all_rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def bytes_per_point(self) -> int:
        """
            Return the memory required per interpolation point *result* in bytes
        """
    @property
    def invA(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def normalizes_vectors(self) -> bool:
        """
            Whether interpolated eigenvectors are normalized, given the stored data; see :py:meth:`set_vector_normalization`
        """
    @property
    def rlu(self) -> numpy.ndarray[numpy.float64]:
        ...
    @property
    def tetrahedra(self) -> list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]]:
        ...
    @property
    def values(self) -> numpy.ndarray[numpy.complex128]:
        """
            Return a shared view of the stored eigenvalues
        """
    @property
    def vector_metric(self) -> list[float]:
        """
            The diagonal metric used to normalize eigenvectors; empty for the identity
        """
    @property
    def vector_normalization(self) -> str:
        """
            When eigenvectors are normalized: ``"automatic"`` (the default), ``"on"`` or ``"off"``
        """
    @property
    def vectors(self) -> numpy.ndarray[numpy.complex128]:
        """
            Return a shared view of the stored eigenvectors
        """
class RotatesLike:
    """
        Enumeration indicating how vector and matrix values transform
      
    
    Members:
    
      vector : Rotates like a vector
    
      pseudovector : Rotates like a pseudovector
    
      Gamma : Rotates like a (real space) phonon eigenvector
    """
    Gamma: typing.ClassVar[RotatesLike]
    __members__: typing.ClassVar[dict[str, RotatesLike]]
    pseudovector: typing.ClassVar[RotatesLike]
    vector: typing.ClassVar[RotatesLike]
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class _LatticeGrid:
    """
        A periodic triangulation of a lattice, invariant under its point group.
    
        Internal: part of the structured mesh under development; for tests only.
        Grid points are integer coordinates in the lattice basis, times `scale`.
      
    """
    scale: typing.ClassVar[int] = 840
    def __init__(self, basis_rows: numpy.ndarray[numpy.float64], operations: list[numpy.ndarray[numpy.int64]]) -> None:
        ...
    def invariant(self, radius: float) -> bool:
        ...
    def locate(self, x: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)]) -> tuple:
        ...
    def patch(self, radius: float) -> numpy.ndarray[numpy.int64]:
        ...
    @property
    def degenerate(self) -> bool:
        ...
    @property
    def pattern(self) -> numpy.ndarray[numpy.int64]:
        ...
class _LatticeBoundary:
    """
        The irreducible zone's boundary, exactly: faces, pairing cells, special points.
    
        Internal: part of the structured mesh under development; for tests only.
        Coordinates are in the primitive reciprocal lattice basis.
      
    """
    def __init__(self, metric: numpy.ndarray[numpy.float64], operations: list[numpy.ndarray[numpy.int64]], cone: list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]] | None = None) -> None:
        ...
    @property
    def cells(self) -> list:
        ...
    @property
    def faces(self) -> list:
        ...
    @property
    def map_count(self) -> int:
        ...
    @property
    def special_points(self) -> numpy.ndarray[numpy.float64]:
        ...
class _LatticeTri:
    """
        The structured mesh of the irreducible zone: the grid clipped to the zone.
    
        Internal: under development; for tests only. Vertices are in the primitive
        reciprocal lattice basis; the grid lattice is that lattice divided by `n`.
      
    """
    def __init__(self, metric: numpy.ndarray[numpy.float64], operations: list[numpy.ndarray[numpy.int64]], n: int, cone: list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]] | None = None) -> None:
        ...
    def planes_of(self, arg0: int) -> list[int]:
        ...
    def refine(self, marked: numpy.ndarray[numpy.int64], min_edge: float = 0.0) -> None:
        """
        Refine the marked tetrahedra (rows of vertex indices) once, with closure
        """
    @property
    def clipped(self) -> int:
        ...
    @property
    def faces(self) -> list:
        ...
    @property
    def self_paired_ties(self) -> int:
        ...
    @property
    def tetrahedra(self) -> numpy.ndarray[numpy.int64]:
        ...
    @property
    def vertices(self) -> numpy.ndarray[numpy.float64]:
        ...
def _lattice_mesh_inputs(brillouin_zone: BrillouinZone) -> tuple:
    """
        (metric, operations, cone, basis) for meshing a zone's irreducible part: the
        primitive reciprocal metric, the point group on primitive reciprocal
        coordinates, the zone's wedge as integer normals (c·x >= 0 inside), and the
        primitive reciprocal vectors as columns. Internal; for tests only.
    """
@typing.overload
def emit() -> bool:
    """
    Return the output status of the :py:mod:`brille` status printer.
    
    Returns
    -------
    bool
        The value of the status printer STDOUT switch.
    """
@typing.overload
def emit(status: bool) -> bool:
    """
    Modify the output status of the :py:mod:`brille` status printer.
    
    Parameters
    ----------
    emt : bool, optional
        Control whether status messages are printed to STDOUT
    
    Returns
    -------
    bool
        The value of the status printer STDOUT switch.
    """
@typing.overload
def emit_datetime() -> bool:
    """
    Return the timestamp output status of the :py:mod:`brille` status printer.
    
    Returns
    -------
    bool
        The value of the status printer timestamp switch
    """
@typing.overload
def emit_datetime(status: bool) -> bool:
    """
    Modify the timestamp output status of the :py:mod:`brille` status printer.
    
    Parameters
    ----------
    emt : bool, optional
        Control whether a timestamp precedes every status message
    
    Returns
    -------
    bool
        The value of the status printer timestamp switch
    """
@typing.overload
def real_space_tolerance() -> float:
    """
    Return the module-global real space floating point tolerance in angstrom.
    """
@typing.overload
def real_space_tolerance(tolerance: float) -> None:
    """
    Set the module-global real space floating point tolerance in angstrom
    """
@typing.overload
def reciprocal_space_tolerance() -> float:
    """
    Return the module-global reciprocal space floating point tolerance in inverse angstrom.
    """
@typing.overload
def reciprocal_space_tolerance(tolerance: float) -> None:
    """
    Set the module-global reciprocal space floating point tolerance in angstrom
    """
__version__: str
build_datetime: str
build_hostname: str
git_branch: str
git_revision: str
version: str
