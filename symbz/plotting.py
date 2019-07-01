from importlib import util
import numpy as np
import matplotlib.pyplot as pp
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import symbz as sbz

def plot_points(x,title=''):
    fig = pp.figure()
    axs = Axes3D(fig)
    axs.scatter(x[:,0],x[:,1],x[:,2],s=10)
    axs.set_title(title)
    pp.show()
def plot_points_with_lines(x,y,title=''):
    fig = pp.figure()
    axs = Axes3D(fig)
    axs.plot(y[:,0],y[:,1],y[:,2])
    axs.scatter(x[:,0],x[:,1],x[:,2],s=10)
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


def plot_bz(bz, origin=None, Q=None, units='invA', face_vectors=False,
            color='b', edgecolor='k', linewidth=1, alpha=0.7):
    if isinstance(bz, (sbz.BZGridQcomplex, sbz.BZGridQ)):
        if Q is None:
            if units == 'rlu':
                Q = bz.mapped_rlu
            elif units == 'invA':
                Q = bz.mapped_invA
        bz = bz.BrillouinZone
    if units == 'rlu':
        verts = bz.vertices
        faces = bz.faces
    elif units == 'primitive':
        verts = bz.vertices_primitive
        faces = bz.faces_primitive
    else:
        verts = bz.vertices_invA
        faces = bz.faces_invA
    if origin is not None and not isinstance(origin, np.ndarray):
        origin = np.array(origin)
    if origin is None or origin.size != 3 or origin.ndim > 1:
        origin = np.array((0, 0, 0))
    # We have the three faces contributing to each vertex
    fpv = bz.faces_per_vertex
    # for plotting we want the opposite -- the vertices contributing per face
    vpf = [(fpv == i).any(1).nonzero()[0] for i in range(len(faces))]
    # if a vertex is on the edge of more than three faces, we miss N-3 here.
    # try to search for any missing vertices for each face:
    for i, f_i in enumerate(faces):
        for j, v_j in enumerate(verts):
            if not (vpf[i] == j).any():
                if np.abs(np.dot(v_j-f_i/2, f_i)) < 1e-15:
                    vpf[i] = np.append(vpf[i], j)

    patches = [np.array([verts[j, :] for j in i]) for i in vpf]
    # patches = [np.array([verts[vpf[i][j], :] for j in range(len(vpf[i]))])
    #            for i in range(len(vpf))]

    # the coordinates of the patch vertices need to be sorted.
    for idx in range(faces.shape[0]):
        # the face vectors are twice the patch-centre vector
        this_patch = patches[idx] - faces[idx]/2
        n_verts = len(this_patch)
        perm = np.array(range(n_verts))
        windings = np.array([
            [vector_angle(x, y, faces[idx]) for x in this_patch]
            for y in this_patch])
        windings[windings == 0] = 1e3
        for i in range(1, n_verts):
            perm[i] = np.argmin(windings[perm[i-1]])
        patches[idx] = patches[idx][perm]+origin

    xyz_min = np.array([x.min() for x in np.vsplit(verts.transpose(), 3)])
    xyz_max = np.array([x.max() for x in np.vsplit(verts.transpose(), 3)])
    dif = xyz_max-xyz_min
    xyz_min -= dif/20
    xyz_max += dif/20

    fig = pp.figure()
    axs = Axes3D(fig)
    collection = Poly3DCollection(patches, edgecolor=edgecolor,
                                  linewidth=linewidth, alpha=alpha)
    collection.set_facecolor(color)
    axs.add_collection3d(collection)

    if face_vectors:
        fvecs = [np.array([[0, 0, 0], faces[i]/2])
                 for i in range(faces.shape[0])]
        lcol = Line3DCollection(fvecs)
        axs.add_collection3d(lcol)
    axs.set_xlim(left=xyz_min[0], right=xyz_max[0])
    axs.set_ylim(bottom=xyz_min[1], top=xyz_max[1])
    axs.set_zlim(bottom=xyz_min[2], top=xyz_max[2])
    if isinstance(Q, np.ndarray) and Q.ndim == 2 and Q.shape[1] == 3:
        axs.scatter(Q[:, 0], Q[:, 1], Q[:, 2])
    pp.show()


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
