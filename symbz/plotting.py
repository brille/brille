from importlib import util
import numpy as np
import matplotlib.pyplot as pp
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

def plot_points(x,title=''):
    fig = pp.figure()
    ax = Axes3D(fig)
    ax.scatter(x[:,0],x[:,1],x[:,2],s=10)
    ax.set_title(title)
    pp.show()
def plot_points_with_lines(x,y,title=''):
    fig = pp.figure()
    ax = Axes3D(fig)
    ax.plot(y[:,0],y[:,1],y[:,2])
    ax.scatter(x[:,0],x[:,1],x[:,2],s=10)
    ax.set_title(title)
    pp.show()

def plot_bz(bz,origin=None,Q=None,units='invA',color='b',edgecolor='k',face_vectors=False,linewidth=1,alpha=0.7):
    if units == 'rlu':
        v = bz.vertices
        f = bz.faces
    elif units == 'primitive':
        v = bz.vertices_primitive
        f = bz.faces_primitive
    else:
        v = bz.vertices_invA
        f = bz.faces_invA
    if origin is not None and type(origin) != np.ndarray:
        origin = np.array(origin)
    if origin is None or origin.size != 3 or origin.ndim>1:
        origin = np.array((0,0,0))
    # We have the three faces contributing to each vertex
    fpv = bz.faces_per_vertex
    # but for plotting we want the opposite -- the vertices contributing per face
    vpf = [ (fpv==i).any(1).nonzero()[0] for i in range(len(f)) ]
    # if a vertex is on the edge of more than three faces, we will miss N-3 here.
    # try to search for any missing vertices for each face:
    for i in range(len(f)):
        for j in range(len(v)):
            if not (vpf[i]==j).any():
                if np.abs(np.dot(v[j]-f[i]/2,f[i])) < 1e-14:
                    vpf[i] = np.append(vpf[i],j)

    patches = [ np.array([v[vpf[i][j],:] for j in range(len(vpf[i]))]) for i in range(len(vpf))]
    # the coordinates of the patch vertices need to be sorted.
    for t in range(f.shape[0]):
        this_patch = patches[t] - f[t]/2 # the face vectors are twice the patch-centre vector
        nv = len(this_patch)
        perm = np.array(range(nv))
        while True:
            tp = this_patch[perm];
            notok = np.array([ np.dot(np.cross(tp[i-1],tp[i]),f[t])<0 for i in range(nv) ])
            if not notok.any():
                break
            n = notok.nonzero()[0]
            if len(n)>1:
                r = n[1]
            else:
                r = np.random.randint(nv)
            perm[[n[0],r]] = perm[[r,n[0]]]
        patches[t] = patches[t][perm]+origin

    min = np.array([x.min() for x in np.vsplit(v.transpose(),3)])
    max = np.array([x.max() for x in np.vsplit(v.transpose(),3)])
    dif = max-min
    min -= dif/20
    max += dif/20

    fig = pp.figure()
    ax = Axes3D(fig)
    collection=Poly3DCollection(patches,edgecolor=edgecolor,linewidth=linewidth,alpha=alpha)
    collection.set_facecolor(color)
    ax.add_collection3d(collection)


    if face_vectors:
        fvecs = [ np.array([[0,0,0],f[i]/2]) for i in range(f.shape[0]) ]
        lcol = Line3DCollection(fvecs)
        ax.add_collection3d(lcol)
    ax.set_xlim(left=min[0],right=max[0])
    ax.set_ylim(bottom=min[1],top=max[1])
    ax.set_zlim(bottom=min[2],top=max[2])
    if Q is not None and type(Q) == numpy.ndarray and Q.ndim==2 and Q.shape[1]==3:
        ax.scatter(Q[:,0],Q[:,1],Q[:,2])
    pp.show()

# def cube(p0,p1):
#     dx = np.array((p1[0]-p0[0],0,0));
#     dy = np.array((0,p1[1]-p0[1],0));
#     dz = np.array((0,0,p1[2]-p0[2]));
#     v = p0+np.array([dx-dx, dx, dx+dy, dy, dz, dz+dx, dz+dx+dy, dz+dy])
#     t = np.array([[0,1,2,3],[0,1,5,4],[0,4,7,3],[4,5,6,7],[6,2,1,5],[2,6,7,3]])
#     patches = [ v[x] for x in t ]
#     return patches
