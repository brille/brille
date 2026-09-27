"""Profile brille's tetrahedral mesh (BZMeshQdc): time, CPU per thread, memory,
mesh quality and interpolation accuracy, stage by stage.

Each (lattice, mesh size, thread count) runs in its own process, so that memory
and threads belong to that case alone. The parent samples the child's threads
from /proc while the child reports its own stage boundaries and measurements.

    python profile/mesh/profile_mesh.py --out results/tetgen
    python profile/mesh/profile_mesh.py --out results/quick --cases fcc --sizes 100 --threads 1

Stages, in order, for one mesh:

bz            BrillouinZone construction (not the mesh; a reference)
build         BZMeshQdc(bz, max_size=V_ir/size, num_levels)
fill          one scalar per vertex: a smooth periodic function with the crystal's symmetry
locate        interpolating that scalar at points already in the irreducible zone:
              dominated by finding each point's tetrahedron
interpolate   the same from points anywhere: adds moving them into the zone
fill_modes    phonon-like data, 12 modes with 12 complex components each
modes         interpolating that data at points already in the irreducible zone
save, load    HDF5 round trip

Per stage the child records wall time, process CPU time, each thread's CPU time
(from /proc/self/task), the resident memory before and after, and the peak
resident memory within the stage (VmHWM, reset through /proc/self/clear_refs).
The parent's samples give how many threads were running over time.

Per mesh it also records vertex and tetrahedron counts, shape quality (minimum
dihedral angle and normalized radius ratio), and the interpolation error of the
scalar function against its exact values.

Linux only.
"""
import argparse
import json
import os
import platform
import selectors
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

#: name: ((lengths, angles), Hall symbol) -- the lattices of wrap/tests/test_20
CASES = {
    "triclinic": (((4.02, 4.90, 3.29), (98.97, 86.24, 88.47)), "-P 1"),
    "monoclinic C": (((5.1, 8.8, 5.2), (90, 104.5, 90)), "-C 2y"),
    "tetragonal I": (((4.0, 4.0, 9.0), (90, 90, 90)), "-I 4 2"),
    "trigonal P 3": (((3.1, 3.1, 4.7), (90, 90, 120)), "P 3"),
    "hexagonal 6/m": (((3.0, 3.0, 5.0), (90, 90, 120)), "-P 6"),
    "rhombohedral -3m": (((4.9, 4.9, 13.8), (90, 90, 120)), '-R 3 2"'),
    "fcc": (((5.64,) * 3, (90,) * 3), "-F 4 2 3"),
}
TICK = os.sysconf("SC_CLK_TCK")
MODES = 12
REAL_LATTICE = 3      # brille's LengthUnit.real_lattice


# --- /proc readers ---------------------------------------------------------------
def thread_times(pid="self"):
    """{thread id: (state, CPU seconds)} for every thread of the process"""
    out = {}
    base = f"/proc/{pid}/task"
    for tid in os.listdir(base):
        try:
            with open(f"{base}/{tid}/stat") as f:
                stat = f.read()
        except OSError:
            continue
        fields = stat[stat.rindex(")") + 2:].split()
        out[int(tid)] = (fields[0], (int(fields[11]) + int(fields[12])) / TICK)
    return out


def memory(pid="self"):
    """(resident, peak resident) in MiB"""
    values = {}
    with open(f"/proc/{pid}/status") as f:
        for line in f:
            if line.startswith(("VmRSS", "VmHWM")):
                key, kb, _ = line.split()
                values[key[:-1]] = int(kb) / 1024
    return values.get("VmRSS", 0.0), values.get("VmHWM", 0.0)


def reset_peak():
    try:
        with open("/proc/self/clear_refs", "w") as f:
            f.write("5")
        return True
    except OSError:
        return False


# --- the child: one mesh, all stages -------------------------------------------------
MARK = "@@profile "


def emit(**record):
    print(MARK + json.dumps(record), flush=True)


class Stage:
    """Measure one stage and report it to the parent"""
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        emit(event="begin", stage=self.name, t=time.monotonic())
        self.peak_reset = reset_peak()
        self.rss0, _ = memory()
        self.threads0 = thread_times()
        self.cpu0 = os.times()
        self.wall0 = time.perf_counter()
        self.extra = {}
        return self

    def __exit__(self, kind, value, traceback):
        wall = time.perf_counter() - self.wall0
        cpu1 = os.times()
        threads1 = thread_times()
        rss1, peak = memory()
        per_thread = sorted((threads1[t][1] - self.threads0.get(t, ("", 0.0))[1] for t in threads1), reverse=True)
        cpu = (cpu1.user - self.cpu0.user) + (cpu1.system - self.cpu0.system)
        emit(event="end", stage=self.name, t=time.monotonic(), ok=kind is None,
             error=None if kind is None else f"{kind.__name__}: {value}",
             wall_s=wall, cpu_s=cpu, parallelism=cpu / wall if wall > 0 else 0.0,
             thread_cpu_s=[round(x, 3) for x in per_thread if x > 0],
             threads=len(threads1), rss0_mib=self.rss0, rss1_mib=rss1,
             peak_mib=peak if self.peak_reset else None, **self.extra)
        return False


def symmetric_function(lattice):
    """A smooth function of q (rlu) with the lattice's periodicity and point symmetry:
    a sum of cos(2π r·q) over the orbits of a few short direct-lattice vectors r"""
    W = np.asarray(lattice.pointgroup.W)
    seeds = [np.array(r) for r in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1), (2, 0, 0))]
    A = np.asarray(lattice.real_vectors, dtype=float).T
    terms = {}
    for r in seeds:
        weight = np.exp(-np.linalg.norm(A.T @ r) ** 2 / 50.0)
        for w in W:
            for s in (1, -1):
                terms[tuple(int(x) for x in s * (w @ r))] = weight
    R = np.array(list(terms), dtype=float)
    weights = np.array(list(terms.values()))

    def f(q):
        return np.cos(2 * np.pi * np.asarray(q) @ R.T) @ weights
    return f


def mesh_quality(vertices, tetrahedra):
    """Minimum dihedral angle (degrees) and normalized radius ratio 3 r_in / R_circ per tetrahedron"""
    p = vertices[tetrahedra]                                    # (n, 4, 3)
    faces = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    normals = []
    areas = []
    for a, b, c in faces:
        n = np.cross(p[:, b] - p[:, a], p[:, c] - p[:, a])
        areas.append(np.linalg.norm(n, axis=1) / 2)
        normals.append(n / np.linalg.norm(n, axis=1)[:, None])
    # the dihedral angle at the edge shared by faces i and j is π minus the angle between outward normals
    centroid = p.mean(axis=1)
    outward = []
    for k, (a, _, _) in enumerate(faces):
        n = normals[k]
        sign = np.sign(np.einsum("ij,ij->i", n, p[:, a] - centroid))
        outward.append(n * sign[:, None])
    dihedral = []
    for i in range(4):
        for j in range(i + 1, 4):
            dihedral.append(np.degrees(np.pi - np.arccos(np.clip(np.einsum("ij,ij->i", outward[i], outward[j]), -1, 1))))
    min_dihedral = np.min(dihedral, axis=0)
    volume = np.abs(np.einsum("ij,ij->i", p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))) / 6
    r_in = 3 * volume / np.sum(areas, axis=0)
    # circumradius from the edge lengths: 24 V R = sqrt((aA+bB+cC)(aA+bB-cC)(aA-bB+cC)(-aA+bB+cC))
    d = lambda i, j: np.linalg.norm(p[:, i] - p[:, j], axis=1)
    aA, bB, cC = d(0, 1) * d(2, 3), d(0, 2) * d(1, 3), d(0, 3) * d(1, 2)
    prod = (aA + bB + cC) * (aA + bB - cC) * (aA - bB + cC) * (-aA + bB + cC)
    R = np.sqrt(np.clip(prod, 0, None)) / (24 * volume)
    ratio = 3 * r_in / R
    q = lambda x: [float(v) for v in np.percentile(x, (0, 1, 5, 50))]
    return dict(min_dihedral_deg_p0_p1_p5_p50=q(min_dihedral), radius_ratio_p0_p1_p5_p50=q(ratio),
                slivers_below_5deg=int(np.sum(min_dihedral < 5)), slivers_below_10deg=int(np.sum(min_dihedral < 10)),
                volume_invA3=float(volume.sum()))


def child(spec):
    import brille

    name, size, threads = spec["case"], spec["size"], spec["threads"]
    rng = np.random.default_rng(spec["seed"])
    lengths_angles, hall = CASES[name]
    lattice = brille.Lattice(lengths_angles, spacegroup=hall)
    f = symmetric_function(lattice)
    n_points = int(spec["points"])
    q = rng.uniform(-1, 1, size=(n_points, 3))
    par = dict(useparallel=threads != 1, threads=threads)
    info = {}

    with Stage("bz"):
        bz = brille.BrillouinZone(lattice)
    ir_volume = bz.ir_polyhedron.volume
    q_ir = np.ascontiguousarray(np.asarray(bz.ir_moveinto(q, threads=max(threads, 0))[0]))
    # the function must be invariant under what brille does to move points
    if not np.allclose(f(q_ir), f(q), atol=1e-9 * np.abs(f(q)).max()):
        emit(event="warning", message="the test function is not invariant under ir_moveinto")

    with Stage("build") as s:
        mesh = brille.BZMeshQdc(bz, max_size=ir_volume / size, num_levels=spec["levels"])
    vertices = np.asarray(mesh.invA)
    tetrahedra = np.asarray(mesh.tetrahedra)
    info.update(vertices=len(vertices), tetrahedra=len(tetrahedra), ir_volume_invA3=ir_volume)
    info.update(mesh_quality(vertices, tetrahedra))
    rlu = np.asarray(mesh.rlu)
    nv = len(rlu)
    scalar = [1, 0, 0, 0, REAL_LATTICE, 0, 0]
    with Stage("fill"):
        mesh.fill(np.ascontiguousarray(f(rlu)[:, None]), np.array(scalar), np.array([1.0, 0.0, 0.0]),
                  np.zeros((nv, 1), dtype=complex), np.array(scalar), np.array([1.0, 0.0, 0.0]), False)
    exact = f(q)
    with Stage("locate") as s:
        at_ir = np.asarray(mesh.ir_interpolate_at(q_ir, do_not_move_points=True, **par)[0])
        s.extra.update(ns_per_point=1e9 * (time.perf_counter() - s.wall0) / n_points)
    with Stage("interpolate") as s:
        got = np.asarray(mesh.ir_interpolate_at(q, **par)[0])
        s.extra.update(ns_per_point=1e9 * (time.perf_counter() - s.wall0) / n_points)
    err = got.reshape(-1) - exact
    scale = np.std(exact)
    info.update(rms_error=float(np.sqrt(np.mean(err ** 2)) / scale), max_error=float(np.abs(err).max() / scale),
                locate_matches_interpolate=bool(np.allclose(at_ir.reshape(-1), got.reshape(-1))))

    values = rng.random((nv, MODES, 1))
    vectors = rng.random((nv, MODES, MODES)) + 1j * rng.random((nv, MODES, MODES))
    with Stage("fill_modes"):
        mesh.fill(values, np.array(scalar), np.array([1.0, 0.0, 0.0]),
                  vectors, np.array([MODES, 0, 0, 0, REAL_LATTICE, 0, 0]), np.array([0.0, 1.0, 0.0]), False)
    with Stage("modes") as s:
        mesh.ir_interpolate_at(q_ir, do_not_move_points=True, **par)
        s.extra.update(ns_per_point=1e9 * (time.perf_counter() - s.wall0) / n_points)

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "mesh.h5")
        with Stage("save") as s:
            mesh.to_file(path)
        info.update(file_mib=os.path.getsize(path) / 2 ** 20)
        del mesh
        with Stage("load"):
            brille.BZMeshQdc.from_file(path)
    emit(event="info", **info)


# --- the parent: run and sample each child ------------------------------------------
def run_case(spec, rate, timeout):
    env = dict(os.environ, BRILLE_NUM_THREADS=str(spec["threads"] or os.cpu_count()),
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    proc = subprocess.Popen([sys.executable, __file__, "--child", json.dumps(spec)], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, env=env, bufsize=1)
    sel = selectors.DefaultSelector()
    sel.register(proc.stdout, selectors.EVENT_READ)
    records, samples, other = [], [], []
    start = time.monotonic()
    period = 1.0 / rate
    buffer = ""
    while True:
        now = time.monotonic()
        if now - start > timeout:
            proc.kill()
            records.append(dict(event="timeout", after_s=timeout))
            break
        try:
            tt = thread_times(proc.pid)
            rss, _ = memory(proc.pid)
            samples.append((now, rss, sum(1 for s, _ in tt.values() if s == "R"), sum(c for _, c in tt.values())))
        except (FileNotFoundError, ProcessLookupError, ValueError):
            pass
        for key, _ in sel.select(timeout=period):
            chunk = os.read(key.fd, 65536).decode()
            if not chunk:
                sel.unregister(proc.stdout)
                break
            buffer += chunk
            *lines, buffer = buffer.split("\n")
            for line in lines:
                if line.startswith(MARK):
                    records.append(json.loads(line[len(MARK):]))
                elif line.strip():
                    other.append(line)
        if proc.poll() is not None and not sel.get_map():
            break
    proc.wait()
    stderr = proc.stderr.read()
    return records, samples, proc.returncode, "\n".join(other[-20:] + [stderr])


def summarize(records, samples):
    """Per stage: the child's measurements plus running-thread statistics from the samples"""
    stages = {}
    begin = {}
    info = {}
    for r in records:
        if r["event"] == "begin":
            begin[r["stage"]] = r["t"]
        elif r["event"] == "end":
            t0, t1 = begin[r["stage"]], r["t"]
            inside = [s for s in samples if t0 <= s[0] <= t1]
            running = [s[2] for s in inside]
            r = dict(r)
            for k in ("event", "t"):
                r.pop(k)
            if running:
                r.update(running_mean=float(np.mean(running)), running_max=int(max(running)),
                         sampled_rss_max_mib=float(max(s[1] for s in inside)), samples=len(inside))
            stages[r.pop("stage")] = r
        elif r["event"] in ("info", "warning", "timeout"):
            info.update({k: v for k, v in r.items() if k != "event"})
    return stages, info


def environment():
    from importlib.metadata import version
    load = os.getloadavg()
    return dict(date=time.strftime("%Y-%m-%dT%H:%M:%S%z"), brille=version("brille"), python=platform.python_version(),
                cpus=os.cpu_count(), load_average_1_5_15=load, machine=platform.machine())


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--child", help=argparse.SUPPRESS)
    p.add_argument("--out", type=Path, help="directory for the JSON results")
    p.add_argument("--label", default="", help="a name for this run (default: the brille version)")
    p.add_argument("--cases", nargs="*", default=list(CASES), choices=list(CASES))
    p.add_argument("--sizes", type=float, nargs="*", default=[100, 1000, 10000],
                   help="max_size = V_ir / size, so size ≈ the number of tetrahedra requested")
    p.add_argument("--levels", type=int, default=3, help="num_levels")
    p.add_argument("--threads", type=int, nargs="*", default=[1, 0], help="0 means all cores")
    p.add_argument("--points", type=float, default=2e5)
    p.add_argument("--rate", type=float, default=50, help="thread samples per second")
    p.add_argument("--timeout", type=float, default=1800, help="seconds per case")
    p.add_argument("--seed", type=int, default=20260927)
    args = p.parse_args(argv)
    if args.child:
        return child(json.loads(args.child))
    if args.out is None:
        p.error("--out is required")
    args.out.mkdir(parents=True, exist_ok=True)
    env = environment()
    label = args.label or env["brille"]
    results = []
    out_file = args.out / f"{label}.json"
    print(f"# {label}: {env}", flush=True)
    for name in args.cases:
        for size in args.sizes:
            for threads in args.threads:
                spec = dict(case=name, size=size, threads=threads, levels=args.levels, points=args.points, seed=args.seed)
                records, samples, code, stderr = run_case(spec, args.rate, args.timeout)
                stages, info = summarize(records, samples)
                results.append(dict(spec=spec, returncode=code, stderr=stderr[-4000:], stages=stages, info=info,
                                    timeline=[(round(t - samples[0][0], 3), round(r, 1), n, round(c, 2)) for t, r, n, c in samples]
                                    if samples else []))
                out_file.write_text(json.dumps(dict(label=label, environment=env, args=vars(args) | {"out": str(args.out)},
                                                    results=results), indent=1, default=str))
                b = stages.get("build", {})
                loc = stages.get("locate", {})
                print(f"{name:18s} 1/{size:<6g} threads {threads or 'all':>3}: "
                      f"{info.get('vertices', '?'):>7} vertices, build {b.get('wall_s', float('nan')):7.2f} s "
                      f"(x{b.get('parallelism', float('nan')):4.1f}, peak {b.get('peak_mib') or float('nan'):7.1f} MiB), "
                      f"locate {loc.get('ns_per_point', float('nan')):8.0f} ns/pt (x{loc.get('parallelism', float('nan')):4.1f}), "
                      f"rms err {info.get('rms_error', float('nan')):.2e}, exit {code}", flush=True)
    env["load_average_after"] = os.getloadavg()
    out_file.write_text(json.dumps(dict(label=label, environment=env, args=vars(args) | {"out": str(args.out)},
                                        results=results), indent=1, default=str))


if __name__ == "__main__":
    sys.exit(main())
