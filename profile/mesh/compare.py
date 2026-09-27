"""Compare profile_mesh.py results: one table per measure, before and after.

    python profile/mesh/compare.py results/tetgen/tetgen.json results/lattice/lattice.json [--markdown out.md]

Cases are matched by lattice, size and thread count. Several result files may be
given for either side (e.g. runs at different sizes) by joining them with commas.
"""
import argparse
import json
from pathlib import Path


def load(paths):
    runs = {}
    label = None
    for path in paths.split(","):
        data = json.loads(Path(path).read_text())
        label = label or data["label"]
        for r in data["results"]:
            s = r["spec"]
            runs[(s["case"], s["size"], s["threads"])] = r
    return label, runs


def stage(r, name, key):
    return r.get("stages", {}).get(name, {}).get(key)


MEASURES = [
    ("vertices", lambda r: r["info"].get("vertices"), "{:.0f}"),
    ("tetrahedra", lambda r: r["info"].get("tetrahedra"), "{:.0f}"),
    ("build (s)", lambda r: stage(r, "build", "wall_s"), "{:.3f}"),
    ("build CPU / wall", lambda r: stage(r, "build", "parallelism"), "{:.1f}"),
    ("build peak RSS (MiB)", lambda r: stage(r, "build", "peak_mib"), "{:.0f}"),
    ("mesh RSS held (MiB)", lambda r: None if stage(r, "build", "rss1_mib") is None
        else stage(r, "build", "rss1_mib") - stage(r, "build", "rss0_mib"), "{:.1f}"),
    ("locate (ns/pt)", lambda r: stage(r, "locate", "ns_per_point"), "{:.0f}"),
    ("locate CPU / wall", lambda r: stage(r, "locate", "parallelism"), "{:.1f}"),
    ("interpolate (ns/pt)", lambda r: stage(r, "interpolate", "ns_per_point"), "{:.0f}"),
    ("modes (ns/pt)", lambda r: stage(r, "modes", "ns_per_point"), "{:.0f}"),
    ("save (s)", lambda r: stage(r, "save", "wall_s"), "{:.3f}"),
    ("load (s)", lambda r: stage(r, "load", "wall_s"), "{:.3f}"),
    ("file (MiB)", lambda r: r["info"].get("file_mib"), "{:.2f}"),
    ("min dihedral (deg)", lambda r: (r["info"].get("min_dihedral_deg_p0_p1_p5_p50") or [None])[0], "{:.2f}"),
    ("slivers < 10 deg", lambda r: r["info"].get("slivers_below_10deg"), "{:.0f}"),
    ("rms error", lambda r: r["info"].get("rms_error"), "{:.2e}"),
    ("rms error × vertices^(2/3)", lambda r: None if r["info"].get("rms_error") is None
        else r["info"]["rms_error"] * r["info"]["vertices"] ** (2 / 3), "{:.2f}"),
]


def fmt(value, spec):
    return "–" if value is None else spec.format(value)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("before")
    p.add_argument("after")
    p.add_argument("--markdown", type=Path)
    args = p.parse_args(argv)
    label_a, a = load(args.before)
    label_b, b = load(args.after)
    keys = [k for k in a if k in b]
    lines = [f"# {label_a} → {label_b}", ""]
    for name, get, spec in MEASURES:
        lines += [f"## {name}", "", f"| lattice | size | threads | {label_a} | {label_b} | ratio |", "|---|---:|---:|---:|---:|---:|"]
        for k in keys:
            x, y = get(a[k]), get(b[k])
            ratio = f"{y / x:.2f}" if isinstance(x, (int, float)) and isinstance(y, (int, float)) and x else "–"
            lines.append(f"| {k[0]} | {k[1]:g} | {k[2] or 'all'} | {fmt(x, spec)} | {fmt(y, spec)} | {ratio} |")
        lines.append("")
    text = "\n".join(lines)
    if args.markdown:
        args.markdown.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
