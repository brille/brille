"""Switch an interpolation from BZTrellisQ to BZMeshQ, and compare the two.

Run ``python trellis_to_mesh.py`` to build both grids for an example lattice at
matching sizes and print their build time, interpolation time and error.  Edit
``lattice`` and ``model`` to compare them for your own system.

The sections between ``# [name]`` and ``# [/name]`` markers are shown in
docs/howto/trellis_to_mesh.rst.
"""
import time

import numpy as np
from brille import BrillouinZone, BZMeshQdc, BZTrellisQdc, Lattice

lattice = Lattice(((4.02, 4.90, 3.29), (98.97, 86.24, 88.47)), spacegroup="-P 1")
bz = BrillouinZone(lattice)


def model(q):
    """A stand-in for your model: one smooth mode with the lattice's periodicity"""
    q = np.atleast_2d(q)
    return np.sum(np.cos(2 * np.pi * q), axis=1)


def grid_data(q):
    values = model(q)[:, np.newaxis]
    vectors = np.zeros((len(values), 1), dtype=complex)
    return values, (1,), vectors, (1,)


points = 5000   # roughly how many vertices each grid should have

# [before]
volume = bz.ir_polyhedron.volume
trellis = BZTrellisQdc(bz, node_volume_fraction=volume / points)
trellis.fill(*grid_data(trellis.rlu))
# [/before]

# [after]
mesh = BZMeshQdc(bz, max_size=volume / (6 * points))
mesh.fill(*grid_data(mesh.rlu))
# [/after]


# [compare]
def measure(make, q, repeat=3):
    """Build a grid, fill it with the model, and time interpolating at q"""
    start = time.perf_counter()
    grid = make()
    built = time.perf_counter() - start
    grid.fill(*grid_data(grid.rlu))
    best = np.inf
    for _ in range(repeat):
        start = time.perf_counter()
        values = grid.ir_interpolate_at(q, threads=1)[0][:, 0]
        best = min(best, time.perf_counter() - start)
    error = np.sqrt(np.mean((values - model(q)) ** 2))
    return len(grid.rlu), built, 1e6 * best / len(q), error


def compare(sizes=(1000, 5000, 20000)):
    """Build a trellis and a mesh of about each size in sizes, and compare them"""
    q = np.random.default_rng(1).uniform(-1, 1, (20000, 3))
    print(f"{'grid':8s} {'vertices':>8s} {'build/s':>8s} {'µs/point':>9s} {'rms error':>10s}")
    for points in sizes:
        for name, make in (("trellis", lambda: BZTrellisQdc(bz, node_volume_fraction=volume / points)),
                           ("mesh", lambda: BZMeshQdc(bz, max_size=volume / (6 * points)))):
            vertices, built, per_point, error = measure(make, q)
            print(f"{name:8s} {vertices:8d} {built:8.2f} {per_point:9.2f} {error:10.2e}")


if __name__ == "__main__":
    compare()
# [/compare]
