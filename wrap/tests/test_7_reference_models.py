#!/usr/bin/env python3
"""Compare brille against reference models evaluated exactly at every q.

The reference models check the physics-critical paths: symmetry unfolding of
eigenvectors, phase conventions, and meshing real lattices. Known failures are
marked ``xfail(strict=True)`` with a reason naming the problem, so a fix shows
up as an unexpected pass.

Crystals sit on real lattices from ``aflow_lattices.json`` and are given to
brille two ways (see :py:func:`harness.brille_lattice`):

* ``hall``: conventional lengths, angles and Hall symbol, with conventional-cell
  eigenvectors;
* ``explicit``: a primitive cell with explicit symmetry operations, as from a
  phonon code.

By default a representative set of space groups is tested; set
``BRILLE_HARNESS_FULL=1`` to sweep all 230 and to mesh every AFLOW lattice.
Run the full sweep with ``BRILLE_NUM_THREADS=1 pytest -n 12`` (pytest-xdist):
each worker process otherwise starts one brille thread per core.
"""
import faulthandler
import json
import os
import subprocess
import sys
import textwrap
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from harness import (
    BornVonKarman,
    aflow_crystal,
    aflow_lattices,
    build,
    compare_modes,
    interior_vertex_indices,
    interpolate,
    point_group_images,
    primitive,
)
from harness.queries import centring_translations, reciprocal_lattice_shifts

FULL = os.environ.get("BRILLE_HARNESS_FULL", "") not in ("", "0")
# every crystal system and centring type, symmorphic and non-symmorphic groups
REPRESENTATIVE = [1, 2, 14, 19, 62, 70, 88, 122, 148, 167, 176, 194, 206, 227]
SPACE_GROUPS = list(range(1, 231)) if FULL else REPRESENTATIVE
ROUTES = ("hall", "explicit")
SEED = 20260923
POINTS_PER_IR = 100
HERE = Path(__file__).parent

# brille holds the GIL throughout its C++ calls, so a Python-level watchdog such as
# pytest-timeout cannot interrupt a hang in brille; faulthandler's C thread can.
TEST_TIMEOUT = float(os.environ.get("BRILLE_HARNESS_TIMEOUT", 300))


@pytest.fixture(autouse=True)
def _hang_watchdog():
    faulthandler.dump_traceback_later(TEST_TIMEOUT, exit=True)
    yield
    faulthandler.cancel_dump_traceback_later()


# (route, ITA number) pairs, for SEED, on which brille fails before the test can
# check anything.
_DUPLICATE_POINT = pytest.mark.xfail(strict=True, raises=RuntimeError, reason="meshing raises 'Duplicate intersection point' (near-coincident polyhedron points)")
_NO_IR_ZONE = pytest.mark.xfail(strict=True, raises=RuntimeError, reason="no irreducible Brillouin zone found for this centred lattice")
_UNREFINED = pytest.mark.xfail(strict=True, raises=RuntimeError, reason="the mesh ignores max_size on the explicit route")
# A hang cannot be xfailed: the watchdog would kill the worker.
_PSEUDO_CUBIC = pytest.mark.skip(reason="hangs: this R lattice is within 6 ppm of fcc, and a 4e-6 zone edge stalls TetGen refinement")
KNOWN_FAILURES = {
    ("hall", 146): _PSEUDO_CUBIC,
    ("explicit", 5): _DUPLICATE_POINT,
    ("explicit", 146): _PSEUDO_CUBIC,
    **{("explicit", n): _NO_IR_ZONE for n in (22, 197, 199, 203)},
    **{("explicit", n): _UNREFINED for n in (127, 200, 201, 215, 221)},
}


def _unfolding_params():
    return [
        pytest.param(route, number, marks=KNOWN_FAILURES.get((route, number), ()), id=f"{route}-{number}")
        for route in ROUTES
        for number in SPACE_GROUPS
    ]


@lru_cache(maxsize=None)
def crystal(number, route="hall"):
    conventional = aflow_crystal(number, np.random.default_rng(SEED + number))
    return conventional if route == "hall" else primitive(conventional)


@lru_cache(maxsize=None)
def filled(number, route="hall", convention="cell"):
    """A filled mesh with at least one interior vertex; thin wedges need a finer mesh."""
    model = BornVonKarman(crystal(number, route), convention=convention)
    for points_per_ir in (POINTS_PER_IR, 3 * POINTS_PER_IR, 10 * POINTS_PER_IR, 30 * POINTS_PER_IR):
        setup = build(model, points_per_ir=points_per_ir, route=route)
        if len(interior_vertex_indices(setup)):
            return setup
    raise RuntimeError(f"no interior mesh vertices for ITA {number}, even at {points_per_ir} points per IR volume")


# Checking one query costs ~(3N)³: ITA 228's 192-atom cell (576 modes) takes ~2 s per query.
# Budget as many queries as a 144-mode cell could afford in 480.
UNFOLDING_BUDGET = 480 * 144**3


def unfolding_queries(setup, number):
    """Every point-group image, shifted by a random reciprocal lattice vector, of interior vertices.

    Uses a random subset of the interior vertices (at least one) when all of them
    would exceed UNFOLDING_BUDGET, so that every operation is still exercised.
    """
    cr = setup.model.crystal
    rng = np.random.default_rng(number)
    vertices = setup.vertices[interior_vertex_indices(setup)]
    operations = len(np.unique(np.asarray(cr.rotations).reshape(-1, 9), axis=0))
    keep = max(1, UNFOLDING_BUDGET // (operations * (3 * cr.natoms) ** 3))
    if len(vertices) > keep:
        vertices = vertices[np.sort(rng.choice(len(vertices), keep, replace=False))]
    return point_group_images(vertices, cr.rotations, rng, max_shift=2, translations=cr.translations)


@pytest.mark.parametrize("number", SPACE_GROUPS)
def test_reference_model_has_crystal_symmetry(number):
    """The oracle itself: eigenvalues are invariant, and eigenvectors periodic, under the space group."""
    cr = crystal(number)
    model = BornVonKarman(cr)
    rng = np.random.default_rng(number)
    q = rng.random((4, 3)) - 0.5
    reference = model.modes(np.repeat(q, len(np.unique(cr.rotations.reshape(-1, 9), axis=0)), axis=0))
    images = model.modes(point_group_images(q, cr.rotations))
    np.testing.assert_allclose(images.values, reference.values, atol=1e-12)
    g = reciprocal_lattice_shifts(rng, len(q), 3, centring_translations(cr.rotations, cr.translations))
    shifted = model.modes(q + g)
    unshifted = model.modes(q)
    assert compare_modes(unshifted, shifted.values, shifted.vectors).worst()["projector_error"] < 1e-8


@pytest.mark.parametrize("route, number", _unfolding_params())
def test_gamma_unfolding(route, number):
    """brille's Γ-table rotation, atom permutation and phase reproduce the exact eigenvectors."""
    setup = filled(number, route)
    q = unfolding_queries(setup, number)
    values, vectors = interpolate(setup, q)
    worst = compare_modes(setup.model.modes(q), values, vectors).worst()
    assert worst["value_error"] < 1e-10
    assert worst["projector_error"] < 1e-8
    assert worst["norm_error"] < 1e-10


@pytest.mark.xfail(strict=True, reason="brille needs eigenvectors periodic in q (the cell phase convention; see docs/phase_convention.rst)")
def test_gamma_unfolding_atom_phase_convention():
    number = 14
    setup = filled(number, convention="atom")
    q = unfolding_queries(setup, number)
    values, vectors = interpolate(setup, q)
    assert compare_modes(setup.model.modes(q), values, vectors).worst()["projector_error"] < 1e-8


@pytest.mark.xfail(strict=True, raises=RuntimeError, reason="time reversal is added as inversion, which has no atom mapping")
def test_time_reversal_without_inversion():
    number = 19
    cr = crystal(number)
    assert not cr.is_centrosymmetric
    model = BornVonKarman(cr)
    setup = build(model, points_per_ir=POINTS_PER_IR, time_reversal=True)
    vertices = setup.vertices[interior_vertex_indices(setup)]
    q = -point_group_images(vertices, cr.rotations)
    values, vectors = interpolate(setup, q)
    assert compare_modes(model.modes(q), values, vectors).worst()["projector_error"] < 1e-8


@pytest.mark.xfail(strict=True, reason="the mesh does not match itself across equivalent zone faces (see GitHub #114)")
def test_continuity_across_zone_faces():
    # triclinic: the first AFLOW P1 lattice is metrically cubic, and a cube's zone meshes consistently
    setup = filled(2)
    vertices = setup.vertices
    moved = np.asarray(setup.brillouin_zone.ir_moveinto(vertices)[0])
    on_face = vertices[~np.all(np.isclose(moved, vertices), axis=1)]
    inside, _ = interpolate(setup, on_face * (1 - 1e-7))
    outside, _ = interpolate(setup, on_face * (1 + 1e-7))
    assert np.abs(inside - outside).max() < 1e-4


RAISED = 3  # exit status of a subprocess whose brille call raised RuntimeError


class BrilleRaised(Exception):
    """A subprocess's brille call raised RuntimeError, rather than crashing or succeeding."""


def run_isolated(script, *args, timeout=300):
    """Run a script in a subprocess. A crash fails; a RuntimeError raises BrilleRaised."""
    result = subprocess.run([sys.executable, "-c", script, *args], capture_output=True, timeout=timeout)
    stderr = result.stderr.decode()[-2000:]
    if result.returncode == RAISED:
        raise BrilleRaised(stderr)
    assert result.returncode == 0, f"exit {result.returncode}\n{stderr}"


_SURFACE_SCRIPT = textwrap.dedent(
    f"""
    import sys
    sys.path.insert(0, {str(HERE)!r})
    import numpy as np
    from harness import BornVonKarman, aflow_crystal, build, interpolate, primitive
    cr = primitive(aflow_crystal(70, np.random.default_rng(0)))
    setup = build(BornVonKarman(cr), points_per_ir=30, route="explicit")
    try:
        interpolate(setup, setup.vertices)
    except RuntimeError as error:
        print(error, file=sys.stderr)
        sys.exit({RAISED})
    """
)


def test_surface_vertices_do_not_crash():
    """Surface points that round off outside the mesh raise instead of segfaulting."""
    try:
        run_isolated(_SURFACE_SCRIPT)
    except BrilleRaised as error:
        assert "not found in tetrahedral mesh" in str(error)


@pytest.mark.xfail(strict=True, raises=BrilleRaised, reason="surface points round off outside the mesh and are not snapped back inside")
def test_surface_vertices_interpolate():
    run_isolated(_SURFACE_SCRIPT)


# Meshing every AFLOW lattice, the way test_3 builds them, in a subprocess each.
# Indices into aflow_lattices.json of lattices that currently fail.
_DUPLICATE = pytest.mark.xfail(strict=True, raises=BrilleRaised, reason="meshing raises 'Duplicate intersection point' (near-coincident polyhedron points)")
AFLOW_MESH_FAILURES = {
    **{i: _DUPLICATE for i in (9, 14, 152, 166, 179, 195, 303, 313, 344, 359, 371, 381)},
    327: pytest.mark.xfail(strict=True, reason="hangs: this R lattice is within 6 ppm of fcc, and a 4e-6 zone edge stalls TetGen refinement"),
    # near-degenerate monoclinic (beta = 90.04 deg): ~5000 vertices instead of ~100, ~40 s alone
    28: pytest.mark.xfail(strict=False, reason="the mesh over-refines near a pseudo-symmetric cell; slow, times out under load"),
}
AFLOW_MESH_TIMEOUT = 240  # below the hang watchdog, so a slow mesh fails rather than kills the worker


def _aflow_params():
    if not FULL:
        return []
    return [
        pytest.param(i, marks=AFLOW_MESH_FAILURES.get(i, ()), id=f"{i}-hall{entry[0]}")
        for i, entry in enumerate(aflow_lattices())
    ]


@pytest.mark.parametrize("index", _aflow_params())
def test_aflow_lattice_meshes(index):
    """Every real lattice can be meshed without crashing or hanging."""
    code = textwrap.dedent(
        """
        import sys, json, brille
        hall, lengths, angles, symbol = json.loads(sys.argv[1])
        bz = brille.BrillouinZone(brille.Lattice((lengths, angles), spacegroup=symbol))
        try:
            brille.BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 100)
        except RuntimeError as error:
            print(error, file=sys.stderr)
            sys.exit(RAISED)
        """
    ).replace("RAISED", str(RAISED))
    entry = json.dumps(aflow_lattices()[index])
    try:
        run_isolated(code, entry, timeout=AFLOW_MESH_TIMEOUT)
    except subprocess.TimeoutExpired:
        pytest.fail(f"meshing {entry} took more than {AFLOW_MESH_TIMEOUT} s")
