"""Regression tests for the ThreadPool that replaced OpenMP.

Crashes run in a subprocess, so that a regression fails the test instead of
killing the test process.
"""
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest
from brille import BrillouinZone, Lattice


def rhombohedral_zone():
    """R-3 in the hexagonal setting, whose primitive transform P is not symmetric."""
    return BrillouinZone(Lattice(((4.9, 4.9, 13.8), (90, 90, 120)), spacegroup="-R 3"))


@pytest.mark.parametrize("threads", [1, 4])
def test_rhombohedral_moveinto(threads):
    bz = rhombohedral_zone()
    Q = (np.random.default_rng(1).random((500, 3)) - 0.5) * 6
    q, tau = bz.moveinto(Q, threads=threads)
    assert np.allclose(q + tau, Q)
    assert all(bz.isinside(q))


@pytest.mark.parametrize("threads", [1, 4])
def test_rhombohedral_ir_moveinto(threads):
    bz = rhombohedral_zone()
    Q = (np.random.default_rng(2).random((500, 3)) - 0.5) * 6
    q, tau, R, invR = bz.ir_moveinto(Q, threads=threads)
    assert np.allclose(np.einsum("nji,nj->ni", R, q) + tau, Q)
    assert all(bz.isinside(q))


def run(code):
    return subprocess.run([sys.executable, "-c", textwrap.dedent(code)], capture_output=True, timeout=300)


def test_mesh_surface_points_do_not_crash():
    """Mesh vertices on the irreducible zone surface can round off outside every tetrahedron.

    TetTri::locate then followed the not-found sentinel into the layer connections and
    segfaulted in a worker thread. Until such points are snapped back inside, a
    RuntimeError naming the missing points is the expected outcome.
    """
    result = run("""
        import numpy as np
        from brille import BrillouinZone, BZMeshQdc, Lattice
        lattice = Lattice(((4.023643, 4.900927, 3.288319), (98.972989, 86.236629, 88.466529)), spacegroup="-P 1")
        bz = BrillouinZone(lattice)
        mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 30)
        mesh.fill(np.ones(len(mesh.rlu)), (1,), mesh.rlu, (0, 3))
        try:
            mesh.ir_interpolate_at(mesh.rlu, threads=4)
        except RuntimeError as error:
            assert "not found in tetrahedral mesh" in str(error), error
        """)
    assert result.returncode == 0, result.stderr.decode()[-2000:]


def test_exception_in_worker_reaches_python():
    """An exception thrown on a worker thread called std::terminate (SIGABRT) instead of raising.

    This monoclinic lattice (beta = 90.04 deg) trips 'Duplicate intersection point' while
    the mesh layers are connected. Once that is fixed, meshing succeeds, which also passes.
    """
    result = run("""
        from brille import BrillouinZone, BZMeshQdc, Lattice
        bz = BrillouinZone(Lattice(((3.580975428, 3.5819754212, 3.5869753869), (90.0, 90.04, 90.0)), spacegroup="P 1 m 1"))
        try:
            BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 100)
        except RuntimeError as error:
            print(error)
        """)
    assert result.returncode == 0, result.stderr.decode()[-2000:]


def pool_workers(setting):
    """Threads started by the first parallel call, with BRILLE_NUM_THREADS set (or unset for None)."""
    code = textwrap.dedent("""
        import os
        import numpy as np
        from brille import BrillouinZone, Lattice
        bz = BrillouinZone(Lattice(((4.9, 4.9, 13.8), (90, 90, 120)), spacegroup="-R 3"))
        before = len(os.listdir("/proc/self/task"))
        bz.moveinto(np.random.default_rng(0).random((100, 3)))
        print(len(os.listdir("/proc/self/task")) - before)
        """)
    env = {k: v for k, v in os.environ.items() if k != "BRILLE_NUM_THREADS"}
    env["OPENBLAS_NUM_THREADS"] = "1"
    if setting is not None:
        env["BRILLE_NUM_THREADS"] = setting
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, timeout=300, env=env)
    assert result.returncode == 0, result.stderr.decode()[-2000:]
    # the last line: importing brille without matplotlib prints a notice to stdout
    return int(result.stdout.decode().split()[-1])


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="counts threads through /proc")
def test_brille_num_threads_sizes_the_pool():
    """BRILLE_NUM_THREADS sets the pool size when no thread count is given; bad values are ignored."""
    assert pool_workers("3") == 3
    assert pool_workers("not-a-number") == pool_workers(None)
