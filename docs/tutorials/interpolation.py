"""Interpolate with brille: a linear field, then an iron spin-wave dispersion.

Run ``python interpolation.py`` to show the figures with matplotlib, or
``python interpolation.py --save DIRECTORY`` to write them as SVG files; the
figures in the documentation were made that way. ``--no-figures`` runs the
calculation alone, as the tests do.

The sections between ``--8<-- [start:name]`` and ``[end:name]`` markers are shown
in docs/tutorials/interpolation.md.
"""
import argparse
from pathlib import Path

import numpy as np

# --8<-- [start:lattice]
from brille import Lattice

lattice = Lattice(((1, 1, 1), (90, 90, 90)), real_space=False)
print(f"reciprocal cell volume {lattice.volume_star:.3f} Å⁻³")
# --8<-- [end:lattice]

# --8<-- [start:zone]
from brille import BrillouinZone

zone = BrillouinZone(lattice)
print(f"first zone {zone.polyhedron.volume:.3f} Å⁻³, irreducible zone {zone.ir_polyhedron.volume:.3f} Å⁻³")
# --8<-- [end:zone]

# --8<-- [start:grid]
from brille import BZMeshQdd

grid = BZMeshQdd(zone, max_size=zone.ir_polyhedron.volume / 100)
print(f"{len(grid.rlu)} grid points")
# --8<-- [end:grid]


# --8<-- [start:fill]
def phi(q):
    """A linear scalar field, the first Cartesian component of q"""
    return np.asarray(q)[:, 0]


values = phi(grid.invA)[:, np.newaxis]   # one value per grid point
elements = (1,)                          # each point holds one scalar
grid.fill(values, elements, values, elements)
# --8<-- [end:fill]

# --8<-- [start:inside]
q = np.random.default_rng(1).random((10, 3)) - 0.5   # inside the first zone
interpolated, _ = grid.ir_interpolate_at(q)
print(f"inside the zone, exact: {np.allclose(interpolated[:, 0], phi(q))}")
# --8<-- [end:inside]
assert np.allclose(interpolated[:, 0], phi(q))

# --8<-- [start:outside]
far = 20 * q                                          # anywhere in (-10, 10)
interpolated, _ = grid.ir_interpolate_at(far)
print(f"outside the zone, phi(Q - G): {np.allclose(interpolated[:, 0], phi(far - np.round(far)))}")
# --8<-- [end:outside]
assert np.allclose(interpolated[:, 0], phi(far - np.round(far)))

# --8<-- [start:iron]
iron = Lattice(((2.87, 2.87, 2.87), (90, 90, 90)), spacegroup="Im-3m")
iron_zone = BrillouinZone(iron)
iron_grid = BZMeshQdd(iron_zone, max_size=iron_zone.ir_polyhedron.volume / 6000)
print(f"iron: irreducible zone {iron_zone.ir_polyhedron.volume / iron_zone.polyhedron.volume:.4f} "
      f"of the first, {len(iron_grid.rlu)} grid points")
# --8<-- [end:iron]


# --8<-- [start:dispersion]
def omega(q, exchange=16, anisotropy=0.01):
    """Iron's acoustic spin-wave energy, in meV, at q in reciprocal lattice units"""
    return anisotropy + 8 * exchange * (1 - np.prod(np.cos(np.pi * np.asarray(q)), axis=1))


energies = omega(iron_grid.rlu)[:, np.newaxis]
iron_grid.fill(energies, elements, energies, elements)
# --8<-- [end:dispersion]

# --8<-- [start:check]
at_vertices, _ = iron_grid.ir_interpolate_at(iron_grid.rlu)
print(f"the grid returns its own values: {np.allclose(at_vertices, energies)}")
# --8<-- [end:check]
assert np.allclose(at_vertices, energies)

# --8<-- [start:path]
corners = np.array([[0, 0, 0], [1, 0, 0], [0.5, 0.5, 0], [0, 0, 0], [0.5, 0.5, 0.5]])
per_leg = 100
path = np.vstack([np.linspace(corners[i], corners[i + 1], per_leg) for i in range(len(corners) - 1)])
along, _ = iron_grid.ir_interpolate_at(path)
error = np.abs(along[:, 0] - omega(path))
print(f"along the path: largest error {error.max():.2f} meV of {omega(path).max():.0f} meV")
# --8<-- [end:path]
assert error.max() < 0.05 * omega(path).max()


def figures(save=None):
    """Plot the interpolated dispersion along the path; save the figures in save, or show them"""
    # --8<-- [start:plot]
    import matplotlib.pyplot as plt

    x = np.arange(len(path))
    ticks = [*(per_leg * np.arange(len(corners) - 1)), len(path) - 1]
    labels = ["(" + " ".join(f"{c:g}" for c in corner) + ")" for corner in corners]

    def plot(ax):
        ax.plot(x, along[:, 0], "-k", label="interpolated")
        ax.plot(x, omega(path), "--r", label="exact")
        ax.set_xticks(ticks, labels)
        ax.set_xlabel(r"$\mathbf{Q}$ (r.l.u.)")
        ax.set_ylabel(r"$\omega(\mathbf{Q})$ (meV)")

    fig, ax = plt.subplots(figsize=(6.4, 3.6), layout="constrained")
    plot(ax)
    ax.legend()
    # --8<-- [end:plot]
    zoom, zax = plt.subplots(figsize=(6.4, 3.6), layout="constrained")
    plot(zax)
    zax.set_xlim(80, 120)
    zax.set_ylim(247, 257)
    if save is None:
        plt.show()
    else:
        fig.savefig(Path(save) / "interpolation_path.svg")
        zoom.savefig(Path(save) / "interpolation_zoom.svg")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--save", metavar="DIRECTORY", help="write the figures as SVG files here")
    parser.add_argument("--no-figures", action="store_true", help="skip the figures")
    args = parser.parse_args()
    if not args.no_figures:
        figures(args.save)
