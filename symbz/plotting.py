"""Provides functionality for drawing first Brillouin zones in 3D."""
import numpy as np
import matplotlib.pyplot as pp
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import symbz as sbz

def _check_axes(axs=None):
    if axs is None:
        if pp.get_fignums() and isinstance(pp.gca(), Axes3D):
            axs = pp.gca()
        else:
            axs = Axes3D(pp.figure())
    return axs


# pylint: disable=c0103
def plot_points(x, axs=None, title=None, show=True):
    """Plot the N points contained in the (N,3) ndarray x."""
    axs = _check_axes(axs)
    axs.scatter(x[:, 0], x[:, 1], x[:, 2], s=10)
    if title is not None:
        axs.set_title(title)
    if show:
        pp.show()

def plot_points_with_lines(x, y, axs=None, title=None, show=True):
    """Plot the N points contained in the (N,3) ndarray x with lines y.

    The M line segments defined by the (M+1,3) ndarray y are drawn before the
    points in x.
    """
    axs = _check_axes(axs)
    axs.plot(y[:, 0], y[:, 1], y[:, 2])
    axs.scatter(x[:, 0], x[:, 1], x[:, 2], s=10)
    if title is not None:
        axs.set_title(title)
    if show:
        pp.show()

# pylint: disable=r0912,r0913,r0914,r0915
def plot_bz(bz, axs=None, origin=None, Q=None, units='invA', irreducible=False,
            face_vectors=False, show=True,
            color='b', edgecolor='k', linewidth=1, alpha=0.7):
    """Plot a BrillouinZone or BZGridQ[complex] object.

    Draw the faces of a first Brillouin zone with color, edgecolor, linewidth,
    and alpha specified. The plotting units are selectable via the keyword
    `units` with valid values 'rlu', 'invA', or 'primitive' and are relative
    lattice units of the reciprocal space spanning lattice, inverse ångstrom,
    or relative lattice units of the primitive reciprocal space spanning
    lattice, respectively. The face vectors defining the Brillouin zone can be
    drawn as well if the keyword `face_vectors` is set to True.

    If a (N,3) numpy.ndarray is provided via the keyword `Q` it will be treated
    as points in the specified units of reciprocal space.
    If a BZGridQ or BZGridQcomplex object is provided and `Q` is omitted or
    None, then the mapped grid points in 'rlu' or 'invA' will be set to `Q`.
    """
    # pylint: disable=no-member
    axs = _check_axes(axs)
    if isinstance(bz, (sbz.BZGridQcomplex, sbz.BZGridQ, sbz.BZMeshQcomplex, sbz.BZMeshQ)):
        if Q is None:
            if units == 'rlu':
                Q = bz.rlu
            elif units == 'invA':
                Q = bz.invA
        bz = bz.BrillouinZone
    if origin is not None and not isinstance(origin, np.ndarray):
        origin = np.array(origin)
    if origin is None or origin.size != 3 or origin.ndim > 1:
        origin = np.array((0, 0, 0))
    # we always draw the 1st Brillouin zone
    if units == 'rlu':
        verts = bz.vertices
    elif units == 'primitive':
        verts = bz.vertices_primitive
    else:
        verts = bz.vertices_invA
    bzcolor = color if not irreducible else "w"
    bzedgecolor = edgecolor if not irreducible else "0.5"
    bzlinestyle = '-' if not irreducible else '--'
    bzalpha = alpha if not irreducible else 0

    # the 1st Brillouin zone has on-face points equal to half the normals
    polybz, xyz_min, xyz_max = _make_poly_collection(verts,
                                                     bz.vertices_per_face,
                                                     origin=origin,
                                                     color=bzcolor,
                                                     edgecolor=bzedgecolor,
                                                     linestyle=bzlinestyle,
                                                     linewidth=linewidth,
                                                     alpha=bzalpha)
    if irreducible:
        if units == 'rlu':
            ir_verts = bz.ir_vertices
        elif units == 'primitive':
            ir_verts = bz.ir_vertices_primitive
        else:
            ir_verts = bz.ir_vertices_invA
        if ir_verts.size > 0:
            polyir, _, _ = _make_poly_collection(ir_verts,
                                                 bz.ir_vertices_per_face,
                                                 origin=origin,
                                                 color=color,
                                                 edgecolor=edgecolor,
                                                 linestyle='-',
                                                 linewidth=linewidth,
                                                 alpha=alpha)
            axs.add_collection3d(polyir)
    axs.add_collection3d(polybz)
    if face_vectors:
        if units == 'rlu':
            norms = bz.normals
            point = bz.points
        elif units == 'primitive':
            norms = bz.normals_primitive
            point = bz.points_primitive
        else:
            norms = bz.normals_invA
            point = bz.points_invA
        fvecs = [np.array([p, p+n]) for p, n in zip(point, norms)]
        lcol = Line3DCollection(fvecs)
        axs.add_collection3d(lcol)
    axs.set_xlim(left=xyz_min[0], right=xyz_max[0])
    axs.set_ylim(bottom=xyz_min[1], top=xyz_max[1])
    axs.set_zlim(bottom=xyz_min[2], top=xyz_max[2])
    if isinstance(Q, np.ndarray) and Q.ndim == 2 and Q.shape[1] == 3:
        axs.scatter(Q[:, 0], Q[:, 1], Q[:, 2])
    axs.set_aspect('equal', 'box')
    if show:
        pp.show()
    return axs

def _make_poly_collection(verts, vpf, origin=None, color='b', edgecolor='k',
                          linestyle='-', linewidth=1, alpha=0.5):
    # vpf lists the ordered vertices which make up each facet
    # for each facet, pick-out the vertices which define its polygon face
    patches = [np.array([verts[j, :] for j in i]) for i in vpf]
    # if an origin has been provided, add it to the patches
    if origin is not None and origin.ndim == 1 and origin.shape[0] == 3:
        for p in patches:
            p += origin
    # find the extent of the patches
    xyz_min = np.array([x.min() for x in np.vsplit(verts.transpose(), 3)])
    xyz_max = np.array([x.max() for x in np.vsplit(verts.transpose(), 3)])
    # plus some nice-for-plotting padding
    dif = xyz_max-xyz_min
    xyz_min -= dif/20
    xyz_max += dif/20
    # and create the collection of polygons in 3D
    collection = Poly3DCollection(patches, edgecolor=edgecolor,
                                  linestyle=linestyle, linewidth=linewidth,
                                  alpha=alpha)
    # which requires that the face color be set after the fact
    collection.set_facecolor(color)
    return (collection, xyz_min, xyz_max)

def __cube(p_0, p_1):
    """Return the patches of a cube bounded by points p_0 and p_1."""
    d_x = np.array((p_1[0]-p_0[0], 0, 0))
    d_y = np.array((0, p_1[1]-p_0[1], 0))
    d_z = np.array((0, 0, p_1[2]-p_0[2]))
    verts = p_0+np.array([d_x-d_x,      # 0 (000)
                          d_x,          # 1 (100)
                          d_x+d_y,      # 2 (110)
                          d_y,          # 3 (010)
                          d_z,          # 4 (001)
                          d_z+d_x,      # 5 (101)
                          d_z+d_x+d_y,  # 6 (111)
                          d_z+d_y])     # 7 (011)
    idx = np.array([[0, 1, 2, 3],   # (000)-(100)-(110)-(010)
                    [0, 1, 5, 4],   # (000)-(100)-(101)-(001)
                    [0, 4, 7, 3],   # (000)-(001)-(011)-(010)
                    [4, 5, 6, 7],   # (001)-(101)-(111)-(011)
                    [6, 2, 1, 5],   # (111)-(110)-(100)-(101)
                    [2, 6, 7, 3]])  # (110)-(111)-(011)-(010)
    patches = [verts[x] for x in idx]
    return patches

def plot_polyhedron(poly, axs=None, setlims=True, show=True, **kwds):
    """Plot a polyhedron"""
    # pylint: disable=no-member
    axs = _check_axes(axs)
    # the 1st Brillouin zone has on-face points equal to half the normals
    coll, xyz_min, xyz_max = _make_poly_collection(poly.vertices,
                                                   poly.vertices_per_face,
                                                   **kwds)
    axs.add_collection3d(coll)
    if setlims:
        axs.set_xlim(left=xyz_min[0], right=xyz_max[0])
        axs.set_ylim(bottom=xyz_min[1], top=xyz_max[1])
        axs.set_zlim(bottom=xyz_min[2], top=xyz_max[2])
    if show:
        pp.show()
    return axs
