"""brille releases the GIL during long calls, and concurrent calls stay correct."""
import threading
import time

import numpy as np
from brille import BrillouinZone, BZMeshQdc, Lattice

GAMMA, ANGSTROM = 2, 1  # RotatesLike, LengthUnit


def filled_grid():
    lattice = Lattice(((3, 4, 5), (90, 90, 90)), spacegroup="P 1", basis=([[0, 0, 0]], [0]))
    bz = BrillouinZone(lattice)
    grid = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 2000)
    q = np.asarray(grid.rlu)
    rng = np.random.default_rng(1)
    vectors = rng.standard_normal((len(q), 3, 3)) + 1j * rng.standard_normal((len(q), 3, 3))
    grid.fill(rng.random((len(q), 3)), (1,), np.ascontiguousarray(vectors), (0, 3, 0, GAMMA, ANGSTROM))
    return grid


def queries(n, seed):
    return np.random.default_rng(seed).random((n, 3)) - 0.5


def test_python_threads_run_during_a_long_interpolation():
    grid = filled_grid()
    q = queries(30_000, 2)
    worker = threading.Thread(target=lambda: grid.ir_interpolate_at(q, threads=1))
    start = time.perf_counter()
    ticks = [start]
    worker.start()
    while worker.is_alive():
        time.sleep(0.001)  # actually sleeps longer on some systems (about 8 ms on macOS runners)
        ticks.append(time.perf_counter())  # pure-Python work, possible only while brille does not hold the GIL
    elapsed = ticks[-1] - start
    worker.join()
    # With the GIL held the main thread cannot run until brille returns: one gap spans
    # the whole call. Released, the gaps are only as long as the sleep resolution.
    assert elapsed > 0.1
    assert max(b - a for a, b in zip(ticks, ticks[1:])) < 0.25 * elapsed


def test_concurrent_calls_match_serial_results():
    grid = filled_grid()
    q = [queries(2_000, seed) for seed in range(4)]
    expected = [grid.ir_interpolate_at(qi, threads=2) for qi in q]
    results = [None] * len(q)
    errors = []

    def work(i):
        try:
            for _ in range(5):
                results[i] = grid.ir_interpolate_at(q[i], useparallel=True, threads=2)
        except Exception as error:  # noqa: BLE001 - reported below
            errors.append(error)

    threads = [threading.Thread(target=work, args=(i,)) for i in range(len(q))]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not errors
    for (ev, evec), (rv, rvec) in zip(expected, results):
        np.testing.assert_array_equal(np.asarray(rv), np.asarray(ev))
        np.testing.assert_array_equal(np.asarray(rvec), np.asarray(evec))


def test_meshing_in_threads_matches_serial():
    lattice = Lattice(((3, 4, 5), (90, 95, 90)), spacegroup="P 1 2 1", basis=([[0, 0, 0]], [0]))
    bz = BrillouinZone(lattice)
    serial = len(BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 500).rlu)
    counts = []
    threads = [threading.Thread(target=lambda: counts.append(len(BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 500).rlu)))
               for _ in range(3)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert counts == [serial] * 3
