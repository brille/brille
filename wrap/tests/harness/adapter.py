"""Build brille objects from a reference model, and call into brille.

Only this module imports brille, so the models and comparisons stay usable
against other implementations.
"""
from dataclasses import dataclass

import numpy as np

from .models import Modes, ReferenceModel

#: Integer values of brille's ``RotatesLike`` enumeration
ROTATES_LIKE = {"vector": 0, "pseudovector": 1, "gamma": 2}
#: Integer values of brille's ``LengthUnit`` enumeration
LENGTH_UNIT = {"none": 0, "angstrom": 1, "inverse_angstrom": 2, "real_lattice": 3, "reciprocal_lattice": 4}


@dataclass
class BrilleSetup:
    """A filled brille grid and what it was built from."""

    model: ReferenceModel
    lattice: object
    brillouin_zone: object
    grid: object
    length_unit: str
    vertex_modes: Modes
    route: str

    @property
    def vertices(self):
        """Grid vertices in reciprocal lattice units, (V, 3)."""
        return np.asarray(self.grid.rlu)


def _to_unit(vectors, lattice, unit):
    """Cartesian (crystal frame) displacement vectors to brille's storage unit."""
    if unit == "real_lattice":
        return vectors @ np.linalg.inv(lattice)  # u = f·A  ⇒  f = u·A⁻¹
    raise NotImplementedError(f"length unit {unit!r} is not handled by the harness yet")


def _from_unit(vectors, lattice, unit):
    if unit == "real_lattice":
        return vectors @ lattice
    raise NotImplementedError(f"length unit {unit!r} is not handled by the harness yet")


def brille_lattice(crystal, route="hall"):
    """A brille ``Lattice`` for ``crystal``, given the way a user would.

    Parameters
    ----------
    crystal : Crystal
    route : {"hall", "explicit"}
        * ``"hall"``: conventional lengths and angles plus the Hall symbol, the
          route ``test_2`` and ``test_3`` use. brille builds the zone from the
          primitive cell internally, but its phonon symmetry table uses the
          conventional atoms and operations (``Lattice::primitive`` leaves them
          unchanged), so eigenvectors must describe the conventional cell.
        * ``"explicit"``: the cell's own vectors plus explicit symmetry
          operations, as for a primitive cell from a phonon code (``test_5``).
    """
    import brille

    basis = (crystal.positions, crystal.types)
    if route == "hall":
        if crystal.hall_symbol is None:
            raise ValueError("the 'hall' route needs a crystal in its conventional setting")
        return brille.Lattice((list(crystal.lengths), list(crystal.angles)), spacegroup=crystal.hall_symbol, basis=basis)
    if route == "explicit":
        return brille.Lattice(crystal.lattice, symmetry=(crystal.rotations, crystal.translations), basis=basis)
    raise ValueError(f"unknown route {route!r}")


def build(model: ReferenceModel, points_per_ir=30, time_reversal=False, length_unit="real_lattice", sort=False,
          route="hall"):
    """Fill a brille mesh with ``model``'s exact modes at every vertex.

    Parameters
    ----------
    model : ReferenceModel
    points_per_ir : float
        The largest tetrahedron volume is the irreducible volume divided by this.
    time_reversal : bool
        Passed to brille's ``BrillouinZone(time_reversal_symmetry=...)``.
    length_unit : str
        Unit the eigenvectors are stored in; a key of :py:data:`LENGTH_UNIT`.
    sort : bool
        Whether brille should determine mode permutations after filling.
    route : {"hall", "explicit"}
        How the lattice is given to brille; see :py:func:`brille_lattice`.
    """
    import brille

    cr = model.crystal
    lat = brille_lattice(cr, route)
    bz = brille.BrillouinZone(lat, time_reversal_symmetry=time_reversal)
    grid = brille.BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / points_per_ir)
    modes = model.modes(np.asarray(grid.rlu))
    n = cr.natoms
    grid.fill(
        np.ascontiguousarray(modes.values[:, :, None]),
        # scalars only, but LengthUnit "none" throws at query time for this RotatesLike
        np.array([1, 0, 0, 0, LENGTH_UNIT["real_lattice"], 0, 0]),
        np.array([1.0, 0.0, 0.0]),
        np.ascontiguousarray(_to_unit(modes.vectors, cr.lattice, length_unit)),
        np.array([0, 3 * n, 0, ROTATES_LIKE[model.rotates_like], LENGTH_UNIT[length_unit], 0, 0]),
        np.array([0.0, 1.0, 0.0]),
        sort,
    )
    return BrilleSetup(model, lat, bz, grid, length_unit, modes, route)


def interpolate(setup: BrilleSetup, q_rlu, threads=1):
    """Interpolate with brille; returns (values (Q, B), Cartesian vectors (Q, B, N, 3))."""
    q = np.ascontiguousarray(np.atleast_2d(q_rlu), dtype=float)
    values, vectors = setup.grid.ir_interpolate_at(q, useparallel=threads != 1, threads=threads)
    nq = len(q)
    values = np.asarray(values).reshape(nq, -1)
    vectors = _from_unit(np.asarray(vectors), setup.model.crystal.lattice, setup.length_unit)
    return values, vectors
