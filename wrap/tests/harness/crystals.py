"""Crystals with real lattices and a prescribed space group.

Lattices come from ``aflow_lattices.json`` (the extract of AFLOW used by
``test_3_aflow_spacegroups.py``), so every test uses a physically realistic
lattice shape. Atoms are added orbit by orbit, one species per orbit, until
spglib confirms the intended space group:

* candidate positions mix random coordinates with special values (0, 1/4, 1/2,
  …, and coordinates tied together such as x, x, 0), so that high-symmetry
  groups can be realised with a small number of atoms;
* each orbit is a new species, so separate orbits never add symmetry;
* atoms closer than a fraction of the mean spacing are rejected.

A crystal is generated in the conventional setting of its Hall number, which
brille can be given directly by Hall symbol, and can be reduced to a primitive
cell with explicit symmetry operations by :py:func:`primitive`.
"""
import json
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np
import spglib

SYMPREC = 1e-4
AFLOW_FILE = Path(__file__).parent.parent / "aflow_lattices.json"


@dataclass(frozen=True)
class Crystal:
    """A crystal and its space group operations.

    Attributes
    ----------
    lattice : (3, 3) array
        Real-space basis vectors as rows, in angstrom.
    positions : (N, 3) array
        Atom positions in fractional coordinates of ``lattice``.
    types : (N,) int array
        Atom species, contiguous from zero.
    rotations : (S, 3, 3) int array
        Rotation parts of the space group operations, acting on fractional
        real-space coordinates as column vectors.
    translations : (S, 3) array
        Translation parts of the operations, fractional.
    number : int
        International Tables space group number.
    hall_number : int
        spglib Hall number of the setting.
    hall_symbol : str or None
        Set for a crystal in the conventional setting of ``hall_number``, which
        brille can be given by symbol; ``None`` for a primitive reduction.
    """

    lattice: np.ndarray
    positions: np.ndarray
    types: np.ndarray
    rotations: np.ndarray
    translations: np.ndarray
    number: int
    hall_number: int
    hall_symbol: Optional[str] = None

    @property
    def natoms(self) -> int:
        return len(self.types)

    @property
    def lengths(self):
        return np.linalg.norm(self.lattice, axis=1)

    @property
    def angles(self):
        """Basis-vector angles (α, β, γ) in degrees."""
        a, l = self.lattice, self.lengths
        return np.degrees([np.arccos(a[j] @ a[k] / (l[j] * l[k])) for j, k in ((1, 2), (0, 2), (0, 1))])

    @property
    def is_centrosymmetric(self) -> bool:
        return any(np.array_equal(r, -np.eye(3, dtype=int)) for r in self.rotations)

    def __repr__(self):
        setting = self.hall_symbol or "primitive"
        return f"Crystal(ITA {self.number}, {setting!r}, {self.natoms} atoms, {len(self.rotations)} ops)"


@lru_cache(maxsize=1)
def aflow_lattices():
    """The AFLOW extract: a list of ``[hall_number, lengths, angles, symbol]``."""
    with open(AFLOW_FILE) as f:
        return json.load(f)


@lru_cache(maxsize=1)
def aflow_by_ita() -> dict:
    """Map each ITA number to the first AFLOW entry for it."""
    out = {}
    for entry in aflow_lattices():
        out.setdefault(spglib.get_spacegroup_type(entry[0]).number, entry)
    return out


@lru_cache(maxsize=1)
def hall_numbers_by_ita() -> dict:
    """Map each ITA number (1–230) to the first spglib Hall number for it."""
    out = {}
    for hall in range(1, 531):
        out.setdefault(spglib.get_spacegroup_type(hall).number, hall)
    return out


def lattice_vectors(lengths, angles):
    """Row basis vectors with a along x and b in the xy plane; angles in degrees."""
    a, b, c = lengths
    al, be, ga = np.radians(angles)
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    return np.array([[a, 0, 0], [b * np.cos(ga), b * np.sin(ga), 0], [cx, cy, np.sqrt(c * c - cx * cx - cy * cy)]])


def _fallback_lattice(hall_number):
    """Lattice parameters in the style of ``test_2_brillouinzone.py`` for groups AFLOW lacks."""
    t = spglib.get_spacegroup_type(hall_number)
    if 143 <= t.number <= 194:
        if t.choice == "R":
            return [6.0, 6.0, 6.0], [70.0, 70.0, 70.0]
        return [5.0, 5.0, 8.0], [90.0, 90.0, 120.0]
    raise ValueError(f"no fallback lattice for ITA {t.number}")


def _orbit(point, rotations, translations):
    images = (np.einsum("sij,j->si", rotations, point) + translations) % 1.0
    keys = np.round(images * 1e6).astype(np.int64) % 1_000_000  # 0.9999999 and 0 share a key
    _, first = np.unique(keys, axis=0, return_index=True)
    return images[np.sort(first)]


def _wrap(dx):
    return dx - np.round(dx)


def _min_distance(lattice, positions):
    """Smallest interatomic distance, including periodic images."""
    shifts = np.array(np.meshgrid(*[[-1, 0, 1]] * 3, indexing="ij")).reshape(3, -1).T
    dx = _wrap(positions[:, None, :] - positions[None, :, :])
    d = np.linalg.norm((dx[:, :, None, :] + shifts[None, None, :, :]) @ lattice, axis=-1)
    d[d < 1e-9] = np.inf
    return d.min()


_SPECIAL = np.array([0, 1 / 8, 1 / 4, 3 / 8, 1 / 2, 5 / 8, 3 / 4, 7 / 8])


def _candidate(rng):
    """A position whose coordinates are random, special, or tied to each other."""
    x = rng.random(3)
    kind = rng.integers(4)
    if kind == 0:  # general position
        return x
    if kind == 1:  # some coordinates special
        return np.where(rng.random(3) < 0.5, rng.choice(_SPECIAL, 3), x)
    if kind == 2:  # tied coordinates, e.g. (x, x, z), (x, -x, z), (x, x, x)
        p = x.copy()
        p[1] = p[0] * rng.choice([1, -1]) + rng.choice([0, 0.5])
        if rng.random() < 0.3:
            p[2] = p[0]
        return p % 1.0
    return rng.choice(_SPECIAL, 3)  # fully special


def crystal_from_lattice(hall_number, lattice, rng, max_atoms=48, min_spacing=0.3, attempts=4000):
    """Add orbits of atoms to ``lattice`` until spglib finds ``hall_number``'s space group.

    Parameters
    ----------
    hall_number : int
        spglib Hall number; its operations define the conventional setting.
    lattice : (3, 3) array
        Row basis vectors, consistent with the setting.
    rng : numpy.random.Generator
    max_atoms : int
        Upper bound on atoms in the conventional cell.
    min_spacing : float
        Reject atoms closer than this fraction of the mean interatomic spacing.
    attempts : int
        Candidate positions to try before giving up.
    """
    ops = spglib.get_symmetry_from_database(hall_number)
    rotations, translations = ops["rotations"], ops["translations"]
    kind = spglib.get_spacegroup_type(hall_number)
    positions, types = np.zeros((0, 3)), np.zeros(0, dtype=int)
    for _ in range(attempts):
        orbit = _orbit(_candidate(rng), rotations, translations)
        if len(positions) + len(orbit) > max_atoms:
            if len(positions) > max_atoms // 2:  # stuck with extra symmetry: start over
                positions, types = np.zeros((0, 3)), np.zeros(0, dtype=int)
            continue
        trial_positions = np.concatenate([positions, orbit])
        trial_types = np.concatenate([types, np.full(len(orbit), types.max(initial=-1) + 1)])
        spacing = (abs(np.linalg.det(lattice)) / len(trial_types)) ** (1 / 3)
        if len(trial_types) > 1 and _min_distance(lattice, trial_positions) < min_spacing * spacing:
            continue
        positions, types = trial_positions, trial_types
        dataset = spglib.get_symmetry_dataset((lattice, positions, types), symprec=SYMPREC)
        if dataset is not None and dataset.number == kind.number:
            return Crystal(lattice, positions, types, rotations, translations, kind.number, hall_number, kind.hall_symbol)
    raise RuntimeError(f"could not realise Hall {hall_number} (ITA {kind.number}) within {max_atoms} atoms")


def aflow_crystal(number, rng, max_atoms=(48, 96, 192)):
    """A conventional crystal of ITA ``number`` on its first AFLOW lattice.

    ``max_atoms`` are tried in turn: some face-centred cubic groups (e.g. 203,
    209, 210, 219, 226, 228) cannot be realised exactly within 48 atoms.
    """
    entry = aflow_by_ita().get(number)
    if entry is None:
        hall = hall_numbers_by_ita()[number]
        lengths, angles = _fallback_lattice(hall)
    else:
        hall, lengths, angles, _ = entry
    for limit in max_atoms[:-1]:
        try:
            return crystal_from_lattice(hall, lattice_vectors(lengths, angles), rng, max_atoms=limit, attempts=1000)
        except RuntimeError:
            pass
    return crystal_from_lattice(hall, lattice_vectors(lengths, angles), rng, max_atoms=max_atoms[-1])


def primitive(crystal: Crystal) -> Crystal:
    """Reduce ``crystal`` to a primitive cell with operations in that setting."""
    cell = (crystal.lattice, crystal.positions, crystal.types)
    lattice, positions, types = (np.asarray(x) for x in spglib.standardize_cell(
        cell, to_primitive=True, no_idealize=True, symprec=SYMPREC))
    sym = spglib.get_symmetry((lattice, positions, types), symprec=SYMPREC)
    return replace(
        crystal,
        lattice=lattice,
        positions=positions,
        types=types.astype(int),
        rotations=np.asarray(sym["rotations"]),
        translations=np.asarray(sym["translations"]),
        hall_symbol=None,
    )
