"""Interface to VisPy for OpenGL visualisation

See `VisPy <https://vispy.org>`_ for installation and configuration directions
"""
from dataclasses import dataclass
from brille import Polyhedron
from vispy.color import Color
from typing import List

@dataclass
class VisPolyhedron:
    """Wrapper to contain visualisation options for a brille.Polyhedron or brille.LPolyhedron

    Attributes
    ----------
    polyhedron : Union[Polyhedron, LPolyhedron]
        The polyhedron to be visualised. A LPolyhedron will be converted automatically to its Cartesian equivalent.
    face_color : Union[str, vispy.color.Color]
        The face color of the visualised polyhedron, default 'black'
    edge_color : Union[str, vispy.color.Color]
        The edge color of the visualised polyhedron, default 'black'
    fill : bool
        Control drawing of the faces, default True
    outline : bool
        Control drawing of the edges, default True

    Note
    ----
    In the case that a polyhedron with vertices expressed in units of a lattice basis is provided, the oriented lattice
    basis vectors will be used to convert the vertices into an orthonormal frame during object initialisation.
    """
    polyhedron: Polyhedron
    face_color: Color = Color('black')
    edge_color: Color = Color('black')
    fill: bool = True
    outline: bool = True
    opacity: float = 0.2

    def __post_init__(self):
        from brille import LPolyhedron
        if isinstance(self.polyhedron, LPolyhedron):
            self.polyhedron = self.polyhedron.to_Cartesian()
        if isinstance(self.face_color, str):
            self.face_color = Color(self.face_color)
        if isinstance(self.edge_color, str):
            self.edge_color = Color(self.edge_color)

    def box(self):
        """Return the minimum and maximum corners of the bounding box"""
        return self.polyhedron.vertices.min(axis=0), self.polyhedron.vertices.max(axis=0)


def vis_polyhedron(polyhedron: VisPolyhedron, **kwargs):
    """Visualise a single wrapped polyhedron

    Parameters
    ----------
    polyhedron : Union[VisPolyhedron, brille.Polyhedron, brille.LPolyhedron]
        The polyhedron to be visualised. Any object with similar `vertices` and `faces` properties may work as well.
    kwargs
        Optional keyword arguments for the VisPolyhedron constructor. Only used if a Polyhedron or LPolyhedron
        is provided as input.
        face_color : Union[str, vispy.color.Color]
            The face color used for Polyhedron or LPolyhedron object input.
        edge_color : Union[str, vispy.color.Color]
            The edge color used for Polyhedron or LPolyhedron object input.
        fill : bool
            Control drawing the faces of the Polyhedron or LPolyhedron object.
        outline : bool
            Control drawing the edges of the Polyhedron or LPolyhedron object.
    """
    if not isinstance(polyhedron, VisPolyhedron):
        polyhedron = VisPolyhedron(polyhedron, **kwargs)
    vis_polyhedra([polyhedron])


def vis_polyhedra(polyhedra: List[VisPolyhedron], **kwargs):
    """Visualise multiple polyhedra

    Parameters
    ----------
    polyhedra : List[Union[VisPolyhedron, brille.Polyhedron, brille.LPolyhedron]]
        The polyhedron to be visualised. Any object with similar `vertices` and `faces` properties may work as well.
    kwargs
        Optional keyword arguments for the VisPolyhedron constructor. Only used if a Polyhedron or LPolyhedron
        is provided in the polyhedron list.
        face_color : Union[List[Union[str, vispy.color.Color]],Union[str, vispy.color.Color]]
            A single face color used for all Polyhedron and LPolyhedron objects or a list of colors which will be tiled
            to the size of the full polyhedra list.
        edge_color : Union[List[Union[str, vispy.color.Color]],Union[str, vispy.color.Color]]
            A single edge color used for all Polyhedron and LPolyhedron objects or a list of colors which will be tiled
            to the size of the full polyhedra list.
        fill : Union[List[bool], bool]
            A single value to control drawing the faces of all Polyhedron and LPolyhedron objects or a list of boolean
            values which will be tiled to the size of the full polyhedra list.
        outline : Union[List[bool], bool]
            A single value to control drawing the edges of all Polyhedron and LPolyhedron objects or a list of boolean
            values which will be tiled to the size of the full polyhedra list.
        opacity : Union[List[float], float]
            A single value to control face opacity of all Polyhedron and LPolyhedron objects or a list of floating point
            values which will be tiled to the size of the full polyhedra list.
    """
    if any([not isinstance(x, VisPolyhedron) for x in polyhedra]):
        n_poly = len(polyhedra)
        face_colors = make_colours(n_poly, color=kwargs.get('face_color', None))
        edge_colors = make_colours(n_poly, color=kwargs.get('edge_color', None))
        fill = make_list(n_poly, kwargs.get('fill', True))
        outline = make_list(n_poly, kwargs.get('outline', True))
        opacity = make_list(n_poly, kwargs.get('opacity', 0.2))
        polyhedra = [p if isinstance(p, VisPolyhedron) else
                     VisPolyhedron(p, face_color=fc, edge_color=ec, fill=f, outline=o)
                     for p, fc, ec, f, o in zip(polyhedra, face_colors, edge_colors, fill, outline)]
    from vispy import scene, app
    bg = Color('white')
    canvas = scene.SceneCanvas(keys='interactive', bgcolor=bg, size=(800, 600))
    view = canvas.central_widget.add_view()
    view.camera = 'arcball'

    for p in polyhedra:
        if p.fill:
            view.add(vis_polyhedron_to_mesh(p.polyhedron, color=p.face_color, opacity=p.opacity))
        if p.outline:
            [view.add(line) for line in vis_polyhedron_boundary(p.polyhedron, color=p.edge_color)]

    # set the camera range based on properties of the polyhedra...
    from numpy import vstack
    boxes = [p.box() for p in polyhedra]
    box_min = vstack([m for m, _ in boxes]).min(axis=0)
    box_max = vstack([m for _, m in boxes]).max(axis=0)
    view.camera.set_range(x=[box_min[0], box_max[0]], y=[box_min[1], box_max[1]], z=[box_min[2], box_max[2]])

    # finalise the VisPy setup and execute
    canvas.show()
    app.run()


def vis_polyhedron_to_mesh(polyhedron: Polyhedron, color=None, opacity=1.0):
    """Construct as VisPy mesh from a Polyhedron object"""
    from vispy.scene.visuals import Mesh
    from numpy import array
    if color is None:
        color = 'red'
    v = polyhedron.vertices.astype('float32')
    faces = polyhedron.faces
    # ensure the faces are all triangulated:
    faces = [[face[0], face[i], face[(i + 1) % len(face)]] for face in faces for i in range(1, len(face) - 1)]
    mesh = Mesh(v, array(faces, dtype='uint32'), shading='flat', color=color)
    mesh.opacity = opacity
    return mesh


def vis_polyhedron_boundary(polyhedron: Polyhedron, color=None):
    """Construct a list of VisPy Line objects from the edges of a Polyhedron object"""
    from vispy.scene import Line
    v, faces = polyhedron.vertices.astype('float32'), polyhedron.faces
    if color is None:
        color = 'black'
    lines = [Line(v[[*face, face[0]]], color=color, method='gl', width=3, antialias=True) for face in faces]
    return lines


def make_colours(n, color=None):
    """Construct a list of colors for use in displaying Polyhedron objects

    Parameters
    ----------
    n : int
        The number of colors required
    color : Union[List[Union[str, vispy.color.Color]], Union[str, vispy.color.Color]]
        If color is not provided, the list of all colors known to VisPy will be used

    Returns
    -------
    List[vispy.color.Color]
        A length-n list of colors. If the starting value of `color` is less than length-n, it will be tiled to length-n.

    Examples
    --------
    >>> make_colours(7, ['red', 'blue', 'green'])
    [Color('red'), Color('blue'), Color('green'), Color('red'), Color('blue'), Color('green'), Color('red')]

    >>> from vispy.color import Color
    >>> make_colors(4, Color('black'))
    [Color('black'), Color('black'), Color('black'), Color('black')]
    """
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


def make_list(n, val):
    """Construct a 1-D list of values for use in displaying Polyhedron objects

    Parameters
    ----------
    n : int
        The number of colors required
    val : Union[List[T], T]

    Returns
    -------
    List[T]
        A length-n list of values. If the starting value is less than length-n, it will be tiled to length-n.
        If the input list is longer than the requested list length, it will be truncated.

    Examples
    --------
    >>> make_list(10, [True, False, True, True])
    [True, False, True, True, True, False, True, True, True, False]

    The items in the input list should all be the same type. In the case of numeric types, the numpy type promotion
    rules will cast them all to a consistent type before tiling or truncating. So an 'unused' input value may influence
    the type of the returned values.

    >>> make_list(2, [2, 101, 3.14159])
    [2.0, 101.0]
    """
    from collections.abc import Iterable
    if isinstance(val, Iterable):
        val = list(val)
    from numpy import ndarray, array, tile
    if not isinstance(val, ndarray):
        val = array(val)
    if val.size < n:
        val = tile(val, 1 + n//val.size)
    return val[:n]
