def vis_polyhedra(polyhedra, colors=None, outline=True):
    from vispy import app, scene
    if colors is None or len(colors) < len(polyhedra):
        colors = make_colours(len(polyhedra), color=colors)

    canvas = scene.SceneCanvas(keys='interactive', bgcolor='white', size=(800, 600))
    view = canvas.central_widget.add_view()
    view.camera = 'arcball'

    meshes = [vis_polyhedron_to_mesh(p, color=c) for p, c in zip(polyhedra, colors)]
    for mesh in meshes:
        view.add(mesh)

    if outline:
        # poly_lines = [vis_polyhedron_boundary(p, color=c) for p, c in zip(polyhedra, colors)]
        poly_lines = [vis_polyhedron_boundary(p) for p in polyhedra]
        for lines in poly_lines:
            for line in lines:
                view.add(line)

    # set the camera range based on properties of the polyhedra...
    #view.camera.set_range(x=[-3, 3])

    canvas.show()


def vis_polyhedron_to_mesh(polyhedron, color=None):
    from vispy.scene.visuals import Mesh
    from numpy import array
    if color is None:
        color = 'red'

    v = polyhedron.vertices.astype('float32')
    faces = polyhedron.faces

    # ensure the faces are all triangulated:
    faces = [[face[0], face[i], face[(i + 1) % len(face)]] for face in faces for i in range(1, len(face) - 1)]

    mesh = Mesh(v, array(faces, dtype='uint32'), shading='flat', color=color)
    mesh.opacity = 0.2

    return mesh


def vis_polyhedron_boundary(polyhedron, color=None):
    from vispy.scene import Line
    v, faces = polyhedron.vertices.astype('float32'), polyhedron.faces

    if color is None:
        color = 'black'

    lines = [Line(v[[*face, face[0]]], color=color, method='gl', width=3, antialias=True) for face in faces]

    return lines


def make_colours(n, color=None, **kwds):
    if color is None:
        from vispy.color import get_color_names
        color = filter(lambda x: len(x) > 1, get_color_names())

    from collections.abc import Iterable
    if isinstance(color, Iterable):
        color = list(color)
    if isinstance(color, str) or (isinstance(color, (list, tuple)) and len(color) == 3):
        color = [color]

    from numpy import ndarray, array, tile
    if not isinstance(color, ndarray):
        color = array(color)
    if color.shape[0] < n:
        color = tile(color, 1+n//color.shape[0])
    return color[0:n]
