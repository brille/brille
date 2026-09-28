"""Build, fill and refine a BZMeshQ grid, drawing each step with VisPy.

Run ``python tutorial_03.py`` to open each figure in an interactive VisPy window,
or ``python tutorial_03.py --save DIRECTORY`` to write them as PNG files instead;
the figures in the documentation were made that way, with ``--backend egl`` to
draw them without a display.  ``--no-figures`` runs the
calculation alone, as the tests do.

The sections between ``# [name]`` and ``# [/name]`` markers are shown in
docs/tutorials/tutorial_03.rst.
"""
import argparse
import tempfile
from pathlib import Path

import numpy as np

# [lattice]
from brille import BrillouinZone, Lattice

lattice = Lattice(((3.0, 3.0, 5.0), (90, 90, 120)), spacegroup="-P 6")
bz = BrillouinZone(lattice)
print(f"first zone volume      {bz.polyhedron.volume:.4f} Å⁻³")
print(f"irreducible volume     {bz.ir_polyhedron.volume:.4f} Å⁻³")
# [/lattice]


# [model]
def acoustic_branch(lattice, scale=10.0):
    """An acoustic-like dispersion, in meV, with the lattice's periodicity and symmetry.

    Sums 1 - cos(2 pi r.q) over the symmetry-equivalents r of the shortest direct
    lattice vectors, so that the energy rises linearly from zero at Gamma.
    """
    rotations = np.asarray(lattice.pointgroup.W)
    star = {tuple(s * (w @ r)) for r in ((1, 0, 0), (0, 0, 1)) for w in rotations for s in (1, -1)}
    star = np.array(sorted(star), dtype=float)

    def energy(q):
        q = np.atleast_2d(q)
        return scale * np.sqrt(np.sum(1 - np.cos(2 * np.pi * q @ star.T), axis=1))

    return energy


energy = acoustic_branch(lattice)
# [/model]


# [data]
def grid_data(q):
    """The values and vectors to give a grid at the points q (rlu).

    One value per point, the energy; this model has no eigenvectors, but the
    grid still expects a vectors array, so give it a single zero per point.
    """
    values = energy(q)[:, np.newaxis]
    vectors = np.zeros((len(values), 1))
    return values, vectors


ELEMENTS = (1,)   # each point holds one scalar, in both values and vectors
# [/data]


# [build]
from brille import BZMeshQdd

mesh = BZMeshQdd(bz, max_size=bz.ir_polyhedron.volume / 200)
values, vectors = grid_data(mesh.rlu)
mesh.fill(values, ELEMENTS, vectors, ELEMENTS)
print(f"initial mesh           {len(mesh.rlu)} vertices, {len(mesh.tetrahedra)} tetrahedra")
# [/build]
initial = [np.asarray(a).copy() for a in (mesh.rlu, mesh.invA, mesh.tetrahedra)]


# [errors]
def tetrahedron_errors(mesh):
    """How far the interpolated energy is from the model at each tetrahedron's centre"""
    rlu = np.asarray(mesh.rlu)
    centres = rlu[np.asarray(mesh.tetrahedra)].mean(axis=1)
    interpolated = mesh.ir_interpolate_at(centres)[0][:, 0]
    return np.abs(interpolated - energy(centres))


def rms_error(mesh, q):
    """The root-mean-square interpolation error at the points q"""
    interpolated = mesh.ir_interpolate_at(q)[0][:, 0]
    return np.sqrt(np.mean((interpolated - energy(q)) ** 2))


test_points = np.random.default_rng(1).uniform(-1, 1, (20000, 3))
errors = tetrahedron_errors(mesh)
print(f"largest centre error   {errors.max():.3f} meV")
print(f"rms error              {rms_error(mesh, test_points):.4f} meV")
# [/errors]


# [refine]
tolerance = 0.05     # meV
resolution = 0.05    # inverse Angstrom
for step in range(10):
    where = tetrahedron_errors(mesh) > tolerance
    points = mesh.refinement_points(where, resolution=resolution)
    if len(points) == 0:
        break
    mesh.refine(where, *grid_data(points), resolution=resolution)
    print(f"step {step}: {where.sum():5d} tetrahedra marked, {len(points):5d} vertices added")
print(f"refined mesh           {len(mesh.rlu)} vertices, {len(mesh.tetrahedra)} tetrahedra")
print(f"largest centre error   {tetrahedron_errors(mesh).max():.3f} meV")
print(f"rms error              {rms_error(mesh, test_points):.4f} meV")
# [/refine]


# [uniform]
uniform = BZMeshQdd(bz, max_size=bz.ir_polyhedron.volume / 200)
values, vectors = grid_data(uniform.rlu)
uniform.fill(values, ELEMENTS, vectors, ELEMENTS)
while len(uniform.rlu) < len(mesh.rlu):
    uniform.refine(None, *grid_data(uniform.refinement_points()))
print(f"uniformly refined      {len(uniform.rlu)} vertices, rms error {rms_error(uniform, test_points):.4f} meV")
# [/uniform]


# [save]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "acoustic.h5"
    mesh.to_file(str(path))
    loaded = BZMeshQdd.from_file(str(path))
print(f"reloaded               {len(loaded.rlu)} vertices, refinable: {loaded.refinable}")
# [/save]
assert np.array_equal(np.asarray(loaded.rlu), np.asarray(mesh.rlu))


# [draw]
def zone_visuals(polyhedron, color, opacity):
    """VisPy visuals for a Brillouin zone polyhedron: translucent faces and an outline"""
    from brille.vis import vis_polyhedron_boundary, vis_polyhedron_to_mesh

    polyhedron = polyhedron.to_Cartesian()
    faces = vis_polyhedron_to_mesh(polyhedron, color=color, opacity=opacity)
    faces.set_gl_state("translucent", depth_test=False)
    return [faces, *vis_polyhedron_boundary(polyhedron, color="black")]


def surface_visuals(vertices, tetrahedra, colors):
    """VisPy visuals for the surface of a tetrahedral mesh: the triangles on its
    boundary, their vertices coloured, and their edges"""
    from vispy.scene.visuals import Line, Mesh

    corners = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
    triangles = np.concatenate([tetrahedra[:, c] for c in corners])
    # a triangle on the boundary belongs to one tetrahedron, one inside to two
    _, first, count = np.unique(np.sort(triangles, axis=1), axis=0, return_index=True, return_counts=True)
    surface = triangles[first[count == 1]]
    faces = Mesh(vertices, surface, vertex_colors=colors)
    # push the faces back a little, so that their edges are drawn over them
    faces.set_gl_state(depth_test=True, polygon_offset_fill=True, polygon_offset=(1, 1))
    edges = np.concatenate([surface[:, [0, 1]], surface[:, [1, 2]], surface[:, [2, 0]]])
    lines = Line(vertices[edges].reshape(-1, 3), connect="segments", color=(0, 0, 0, 0.5), width=1)
    return [faces, lines]


def gamma_label(size=32):
    """A label for the zone centre"""
    from vispy.scene.visuals import Text

    return [Text("Γ", pos=(0, 0, 0), color="black", font_size=size, anchor_x="right", anchor_y="top")]


def energy_colors(q):
    """Colours for the model energy at q, from dark (zero) to yellow (the maximum)"""
    from vispy.color import get_colormap

    e = energy(q)
    return get_colormap("viridis").map(e / e.max())


def show(visuals, around, save=None, azimuth=-60, elevation=25, size=(640, 520)):
    """Draw visuals in a new VisPy canvas, looking at the points around; show it
    interactively or write it to a PNG file"""
    from vispy import app, io, scene

    canvas = scene.SceneCanvas(keys="interactive", bgcolor="white", size=size, show=save is None)
    view = canvas.central_widget.add_view()
    low, high = np.min(around, axis=0), np.max(around, axis=0)
    view.camera = scene.TurntableCamera(azimuth=azimuth, elevation=elevation, fov=30,
                                        center=tuple((low + high) / 2), scale_factor=np.linalg.norm(high - low))
    for visual in visuals:
        view.add(visual)
    if save is None:
        app.run()
    else:
        io.write_png(str(save), canvas.render())
# [/draw]


def figures(save=None):
    """Draw the tutorial's figures; save them in the directory save, or show them"""
    def target(name):
        return None if save is None else Path(save) / f"tutorial_03_{name}.png"

    # [figure-zone]
    show(zone_visuals(bz.polyhedron, "lightblue", 0.2) + zone_visuals(bz.ir_polyhedron, "orange", 0.6)
         + gamma_label(), bz.polyhedron.to_Cartesian().vertices, save=target("zone"))
    # [/figure-zone]
    # [figure-mesh]
    rlu, xyz, tetrahedra = initial
    show(surface_visuals(xyz, tetrahedra, energy_colors(rlu)) + gamma_label(), xyz, save=target("mesh"))
    # [/figure-mesh]
    # [figure-refined]
    rlu, xyz = np.asarray(mesh.rlu), np.asarray(mesh.invA)
    show(surface_visuals(xyz, np.asarray(mesh.tetrahedra), energy_colors(rlu)) + gamma_label(), xyz, save=target("refined"))
    # [/figure-refined]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--save", metavar="DIRECTORY", help="write the figures as PNG files here")
    parser.add_argument("--backend", help="the VisPy backend, e.g. egl to draw without a display")
    parser.add_argument("--no-figures", action="store_true", help="skip the figures")
    args = parser.parse_args()
    if not args.no_figures:
        if args.backend:
            import vispy
            vispy.use(app=args.backend)
        figures(args.save)
