"""Time each stage of brille's interpolation pipeline separately (DEFICIENCIES.md #27).

Run from ``wrap/tests``::

    python -m harness.benchmark                      # rock salt, a few mesh sizes
    python -m harness.benchmark --points 1e7 --threads 1 0
    python -m harness.benchmark --ita 62 --route explicit --trellis 0.01 0.001

Stages, each in nanoseconds per q point:

move
    ``BrillouinZone.ir_moveinto``: reduce to the first zone, then find the
    point-group operation that lands in the irreducible wedge.
locate
    Interpolating one scalar per vertex at already-moved points. The data is
    trivial, so this is dominated by finding the containing tetrahedron.
interp
    Extra cost of interpolating the full eigenvector data, stored as plain
    complex scalars so that nothing is rotated.
rotate
    Extra cost of storing the same data as Γ-transforming vectors. brille still
    applies the (identity) operation to points that were not moved, so this is
    the rotation cost.
other
    Full ``ir_interpolate_at`` from unmoved points, minus all of the above:
    binding overhead, and anything the stage split missed.
"""
import argparse
import time

import numpy as np
import spglib

from .adapter import LENGTH_UNIT, ROTATES_LIKE, brille_lattice
from .crystals import Crystal, aflow_crystal, primitive
from .models import BornVonKarman


def rock_salt(a=5.64):
    """A two-atom fcc crystal with the full m-3m point group (48 operations)."""
    lattice = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    positions = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    types = np.array([0, 1])
    sym = spglib.get_symmetry((lattice, positions, types))
    return Crystal(lattice, positions, types, sym["rotations"], sym["translations"], 225, 523)


def _best_time(fn, repeat):
    best = np.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def _fill(grid, values, value_elements, vectors, vector_elements):
    grid.fill(
        np.ascontiguousarray(values),
        np.array(value_elements),
        np.array([1.0, 0.0, 0.0]),
        np.ascontiguousarray(vectors),
        np.array(vector_elements),
        np.array([0.0, 1.0, 0.0]),
        False,
    )


def _grids(make_grid, model):
    """Three copies of one mesh: scalar-only, vectors-as-scalars, and Γ vectors."""
    real = LENGTH_UNIT["real_lattice"]
    scalar_elements = [1, 0, 0, 0, real, 0, 0]

    g_scalar = make_grid()
    nv = len(np.asarray(g_scalar.rlu))
    modes = model.modes(np.asarray(g_scalar.rlu))
    nb, n = modes.vectors.shape[1:3]
    values = modes.values[:, :, None]
    vectors = modes.vectors @ np.linalg.inv(model.crystal.lattice)
    _fill(g_scalar, np.ones((nv, 1, 1)), scalar_elements, np.zeros((nv, 1, 1), dtype=complex), scalar_elements)

    g_plain = make_grid()
    _fill(g_plain, values, scalar_elements, vectors.reshape(nv, nb, 3 * n), [3 * n, 0, 0, 0, real, 0, 0])

    g_gamma = make_grid()
    _fill(g_gamma, values, scalar_elements, vectors, [0, 3 * n, 0, ROTATES_LIKE["gamma"], real, 0, 0])
    return g_scalar, g_plain, g_gamma


def run(make_grid, bz, model, q, threads, repeat):
    t_build = time.perf_counter()
    g_scalar, g_plain, g_gamma = _grids(make_grid, model)
    t_build = (time.perf_counter() - t_build) / 3

    par = dict(useparallel=threads != 1, threads=threads)
    move_threads = 1 if threads == 1 else (threads or 0)  # ir_moveinto: 0 means all
    q_ir = np.ascontiguousarray(np.asarray(bz.ir_moveinto(q, threads=move_threads)[0]))
    no_move = dict(do_not_move_points=True, **par)

    t_move = _best_time(lambda: bz.ir_moveinto(q, threads=move_threads), repeat)
    t_scalar = _best_time(lambda: g_scalar.ir_interpolate_at(q_ir, **no_move), repeat)
    t_plain = _best_time(lambda: g_plain.ir_interpolate_at(q_ir, **no_move), repeat)
    t_gamma = _best_time(lambda: g_gamma.ir_interpolate_at(q_ir, **no_move), repeat)
    t_full = _best_time(lambda: g_gamma.ir_interpolate_at(q, **par), repeat)

    per = 1e9 / len(q)
    return {
        "vertices": len(np.asarray(g_gamma.rlu)),
        "build_s": t_build,
        "move": t_move * per,
        "locate": t_scalar * per,
        "interp": (t_plain - t_scalar) * per,
        "rotate": (t_gamma - t_plain) * per,
        "other": (t_full - t_move - t_gamma) * per,
        "full": t_full * per,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--ita", type=int, default=0, help="AFLOW-lattice crystal with this space group (default: rock salt)")
    p.add_argument("--route", choices=("hall", "explicit"), default="hall",
                   help="how an --ita crystal is given to brille (rock salt is always explicit)")
    p.add_argument("--points", type=float, default=1e6, help="number of random q points")
    p.add_argument("--extent", type=float, default=3.0, help="q drawn uniformly from [-extent, extent]³ rlu")
    p.add_argument("--mesh", type=float, nargs="*", default=[100, 1000, 10000],
                   help="BZMesh sizes, as tetrahedra-volume divisors of the irreducible volume")
    p.add_argument("--levels", type=int, default=3, help="BZMesh num_levels")
    p.add_argument("--trellis", type=float, nargs="*", default=[],
                   help="BZTrellis node_volume_fraction values to compare against")
    p.add_argument("--threads", type=int, nargs="*", default=[1, 0], help="thread counts; 0 means all")
    p.add_argument("--repeat", type=int, default=3)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args(argv)

    import brille

    rng = np.random.default_rng(args.seed)
    if args.ita == 0:
        cr, route = rock_salt(), "explicit"
    else:
        cr, route = aflow_crystal(args.ita, rng), args.route
        cr = cr if route == "hall" else primitive(cr)
    model = BornVonKarman(cr)
    lat = brille_lattice(cr, route)
    bz = brille.BrillouinZone(lat)
    ir_volume = bz.ir_polyhedron.volume
    q = rng.uniform(-args.extent, args.extent, size=(int(args.points), 3))

    grids = [(f"mesh 1/{d:g}", lambda d=d: brille.BZMeshQdc(bz, max_size=ir_volume / d, num_levels=args.levels))
             for d in args.mesh]
    grids += [(f"trellis {f:g}", lambda f=f: brille.BZTrellisQdc(bz, node_volume_fraction=f))
              for f in args.trellis]

    print(f"# {cr!r}, {len(q):.0e} points in [-{args.extent:g}, {args.extent:g}]³ rlu, best of {args.repeat}")
    print("| grid | vertices | build (s) | threads | move | locate | interp | rotate | other | full (ns/pt) | 10⁹ pts (s) |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for name, make_grid in grids:
        for threads in args.threads:
            r = run(make_grid, bz, model, q, threads, args.repeat)
            print(f"| {name} | {r['vertices']} | {r['build_s']:.2f} | {threads or 'all'} | "
                  + " | ".join(f"{r[k]:.0f}" for k in ("move", "locate", "interp", "rotate", "other", "full"))
                  + f" | {r['full']:.0f} |", flush=True)


if __name__ == "__main__":
    main()
