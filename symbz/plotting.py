"""Provides functionality for drawing first Brillouin zones in 3D."""
import numpy as np
import matplotlib.pyplot as pp
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import symbz as sbz

# pylint: disable=c0103
def plot_points(x, title=''):
    """Plot the N points contained in the (N,3) ndarray x."""
    fig = pp.figure()
    axs = Axes3D(fig)
    axs.scatter(x[:, 0], x[:, 1], x[:, 2], s=10)
    axs.set_title(title)
    pp.show()

def plot_points_with_lines(x, y, title=''):
    """Plot the N points contained in the (N,3) ndarray x with lines y.

    The M line segments defined by the (M+1,3) ndarray y are drawn before the
    points in x.
    """
    fig = pp.figure()
    axs = Axes3D(fig)
    axs.plot(y[:, 0], y[:, 1], y[:, 2])
    axs.scatter(x[:, 0], x[:, 1], x[:, 2], s=10)
    axs.set_title(title)
    pp.show()

def vector_angle(v_0, v_1, normal=np.array([0, 0, 1])):
    """Calculate the angle between two vectors, given a vector normal."""
    hat_0 = v_0/np.sqrt(np.dot(v_0, v_0))
    hat_1 = v_1/np.sqrt(np.dot(v_1, v_1))
    dot01 = np.dot(hat_0, hat_1)
    crs01 = np.cross(hat_0, hat_1)
    x_vec = dot01*hat_0
    y_vec = hat_1 - x_vec
    y_len = np.sqrt(np.dot(y_vec, y_vec)) * np.sign(np.dot(crs01, normal))
    ang01 = np.arctan2(y_len, dot01)
    return ang01 if ang01 >= 0 else ang01 + 2*np.pi

# pylint: disable=r0912,r0913,r0914,r0915
def plot_bz(bz, origin=None, Q=None, units='invA', irreducible=False,
            face_vectors=False,
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
    fig = pp.figure()
    axs = Axes3D(fig)
    if isinstance(bz, (sbz.BZGridQcomplex, sbz.BZGridQ)):
        if Q is None:
            if units == 'rlu':
                Q = bz.mapped_rlu
            elif units == 'invA':
                Q = bz.mapped_invA
        bz = bz.BrillouinZone
    if origin is not None and not isinstance(origin, np.ndarray):
        origin = np.array(origin)
    if origin is None or origin.size != 3 or origin.ndim > 1:
        origin = np.array((0, 0, 0))
    # we always draw the 1st Brillouin zone
    if units == 'rlu':
        verts = bz.vertices
        norms = bz.faces
    elif units == 'primitive':
        verts = bz.vertices_primitive
        norms = bz.faces_primitive
    else:
        verts = bz.vertices_invA
        norms = bz.faces_invA
    bzcolor = color if not irreducible else "w"
    bzedgecolor = edgecolor if not irreducible else "0.75"
    bzlinestyle = '-' if not irreducible else '--'
    bzalpha = alpha if not irreducible else 0

    # the 1st Brillouin zone has on-face points equal to half the normals
    polybz, xyz_min, xyz_max = _make_poly_collection(verts, norms, norms/2,
                                                     bz.faces_per_vertex,
                                                     origin=origin,
                                                     color=bzcolor,
                                                     edgecolor=bzedgecolor,
                                                     linestyle=bzlinestyle,
                                                     linewidth=linewidth,
                                                     alpha=bzalpha)
    if irreducible:
        if units == 'rlu':
            ir_verts = bz.ir_vertices
            ir_norms = bz.ir_face_normals
            ir_point = bz.ir_face_points
        elif units == 'primitive':
            ir_verts = bz.ir_vertices_primitive
            ir_norms = bz.ir_face_normals_primitive
            ir_point = bz.ir_face_points_primitive
        else:
            ir_verts = bz.ir_vertices_invA
            ir_norms = bz.ir_face_normals_invA
            ir_point = bz.ir_face_points_invA
        if ir_verts.size > 0:
            polyir, _, _ = _make_poly_collection(ir_verts, ir_norms, ir_point,
                                                 bz.ir_faces_per_vertex,
                                                 origin=origin,
                                                 color=color,
                                                 edgecolor=edgecolor,
                                                 linestyle='-',
                                                 linewidth=linewidth,
                                                 alpha=alpha)
            axs.add_collection3d(polyir)
    axs.add_collection3d(polybz)
    if face_vectors:
        fvecs = [np.array([[0, 0, 0], p]) for p in point]
        # fvecs = [np.array([[0, 0, 0], point[i]])
        #          for i in range(point.shape[0])]
        lcol = Line3DCollection(fvecs)
        axs.add_collection3d(lcol)
    axs.set_xlim(left=xyz_min[0], right=xyz_max[0])
    axs.set_ylim(bottom=xyz_min[1], top=xyz_max[1])
    axs.set_zlim(bottom=xyz_min[2], top=xyz_max[2])
    if isinstance(Q, np.ndarray) and Q.ndim == 2 and Q.shape[1] == 3:
        axs.scatter(Q[:, 0], Q[:, 1], Q[:, 2])
    axs.set_aspect('equal','box')
    pp.show()
    return axs

def _make_poly_collection(verts, norms, point, fpv, origin=None,
                          color='b', edgecolor='k',
                          linestyle='-', linewidth=1, alpha=0.5):
    # We have the three faces contributing to each vertex
    # for plotting we want the opposite -- the vertices contributing per face
    vpf = [(fpv == i).any(1).nonzero()[0] for i in range(len(norms))]
    # if a vertex is on the edge of more than three faces, we miss N-3 here.
    # try to search for any missing vertices for each face:
    for i, (n_i, p_i) in enumerate(zip(norms, point)):
        for j, v_j in enumerate(verts):
            if not (vpf[i] == j).any():
                if np.abs(np.dot(v_j-p_i, n_i)) < 1e-15:
                    vpf[i] = np.append(vpf[i], j)

    patches = [np.array([verts[j, :] for j in i]) for i in vpf]
    # the coordinates of the patch vertices need to be sorted.
    for idx in range(norms.shape[0]):
        # find the patch points relative to their centroid
        this_patch = patches[idx] - np.mean(patches[idx], axis=0)
        n_verts = len(this_patch)
        perm = np.array(range(n_verts))
        windings = np.array([
            [vector_angle(x, y, norms[idx]) for x in this_patch]
            for y in this_patch])
        windings[windings == 0] = 1e3
        # or just set the diagonal to a large value?
        # windings[np.diag_indices(n_verts)] = 1e3
        for i in range(1, n_verts):
            perm[i] = np.argmin(windings[perm[i-1]])
        patches[idx] = patches[idx][perm]+origin

    xyz_min = np.array([x.min() for x in np.vsplit(verts.transpose(), 3)])
    xyz_max = np.array([x.max() for x in np.vsplit(verts.transpose(), 3)])
    dif = xyz_max-xyz_min
    xyz_min -= dif/20
    xyz_max += dif/20
    collection = Poly3DCollection(patches, edgecolor=edgecolor,
                                  linestyle=linestyle, linewidth=linewidth,
                                  alpha=alpha)
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
