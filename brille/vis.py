"""Interface to VisPy for OpenGL visualisation

See `VisPy <https://vispy.org>`_ for installation and configuration directions
"""
from dataclasses import dataclass, field
from brille import Polyhedron
from vispy.color import Color
from typing import List

@dataclass
class VisPolyhedron:
    """A wrapper class to contain visualisation options for a :py:class:`brille.Polyhedron` or :py:class:`brille.LPolyhedron`

    Attributes
    ----------
    polyhedron : Union[:py:class:`~brille.Polyhedron`, :py:class:`~brille.LPolyhedron`]
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
    face_color: Color = field(default_factory=lambda: Color('black'))
    edge_color: Color = field(default_factory=lambda: Color('black'))
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
    """Visualise a single :py:class:`~brille.vis.VisPolyhedron` wrapped polyhedron

    Parameters
    ----------
    polyhedron : Union[:py:class:`~brille.vis.VisPolyhedron`, :py:class:`~brille.Polyhedron`, :py:class:`~brille.LPolyhedron`]
        The polyhedron to be visualised. Any object with similar `vertices` and `faces` properties may work as well.
    **kwargs
        Optional keyword arguments for the VisPolyhedron constructor. Only used if a Polyhedron or LPolyhedron
        is provided as input. See other parameters below for expected keyword arguments.

    Other Parameters
    ----------------
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

    This function can be used to visualise the polyhedra that make up a space-filling mesh, but doing so will likely be slow.
    Pre-filtering out all polyhedra which can not be seen (those fully surrounded by other polyhedra) may improve the visualisation speed.


    The visualisation includes interactive elements:
        - The view can be changed via the mouse
        - The shading-style can be cycled by pressing 's' between (`None`, `'smooth'`, and `'flat')
        - The shading `shininess` can be increased by pressing '+' or decreased by pressing '-'
        - The polygon facet fill and wire display can be cycled by pressing 'f'
        - The wireframe display can be toggled by pressing 'w'


    Parameters
    ----------
    polyhedra : List[Union[:py:class:`~brille.vis.VisPolyhedron`, :py:class:`~brille.Polyhedron`, :py:class:`~brille.LPolyhedron`]]
        The polyhedron to be visualised. Any object with similar `vertices` and `faces` properties may work as well.
    **kwargs
        Optional keyword arguments for the VisPolyhedron constructor. Only used if a Polyhedron or LPolyhedron
        is provided in the polyhedron list.

    Other Parameters
    ----------------
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
    elevation: number
        A single value to control the view 'elevation', the angle in degrees above the (x, y) plane of the view
    azimuth: number
        A single value to control the view 'azimuth', the angle in degrees (from -1 * \hat(y)?) of the view
    wireframe_color: Union[str, vispy.color.Color]
        A single color used to draw the wireframe of *every* polyhedron.
    axis: bool
        Control whether a set of axes are displayed along x, y, and z
    """
    if is_interactive():
        print("Visualisation via VisPy may not work properly in an interactive Python shell.")
        print("If no window appears, try embedding your work in a Python script and execute it non-interactively.")

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

    cam = scene.TurntableCamera(elevation=kwargs.get('elevation', 15), azimuth=kwargs.get('azimuth', 150)) 
    view.camera = cam

    if kwargs.get('axis', False):
        axis = scene.visuals.XYZAxis(parent=view.scene)

    wireframe_filter, shading_filter, attach_headlight = make_filters(0, 1, kwargs.get('wireframe_color', 'black'))
    for p in polyhedra:
        if p.fill:
            mesh = vis_polyhedron_to_mesh(p.polyhedron, color=p.face_color, opacity=p.opacity, shading=shading_filter, wireframe=wireframe_filter)
            view.add(mesh)
            add_shading_wireframe_interface(canvas, mesh, shading_filter, wireframe_filter)
        if p.outline:
            [view.add(line) for line in vis_polyhedron_boundary(p.polyhedron, color=p.edge_color)]
    
    attach_headlight(view)

    # set the camera range based on properties of the polyhedra...
    from numpy import vstack
    boxes = [p.box() for p in polyhedra]
    box_min = vstack([m for m, _ in boxes]).min(axis=0)
    box_max = vstack([m for _, m in boxes]).max(axis=0)
    view.camera.set_range(x=[box_min[0], box_max[0]], y=[box_min[1], box_max[1]], z=[box_min[2], box_max[2]])

    # finalise the VisPy setup and execute
    canvas.show()
    app.run()


def vis_polyhedron_to_mesh(polyhedron: Polyhedron, color=None, opacity=1.0, shading=None, wireframe=None):
    """Construct as VisPy mesh from a Polyhedron object"""
    from vispy.scene.visuals import Mesh
    from numpy import array
    if color is None:
        color = 'red'
    v = polyhedron.vertices.astype('float32')
    faces = polyhedron.faces
    # ensure the faces are all triangulated:
    faces = [[face[0], face[i], face[(i + 1) % len(face)]] for face in faces for i in range(1, len(face) - 1)]
    mesh = Mesh(v, array(faces, dtype='uint32'), color=color)
    mesh.opacity = opacity
    if wireframe is not None:
        mesh.attach(wireframe)
    if shading is not None:
        mesh.attach(shading)
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


def make_filters(shininess, width, wireframe_color):
    from vispy.visuals.filters import ShadingFilter, WireframeFilter
    shading_filter = ShadingFilter()
    shading_filter.shading=None
    wireframe_filter = WireframeFilter(width=width, color=wireframe_color)

    def attach_headlight(view):
        light_dir = (0, 1, 0, 0)
        shading_filter.light_dir = light_dir[:3]
        initial_light_dir = view.camera.transform.imap(light_dir)

        @view.scene.transform.changed.connect
        def on_transform_change(event):
            transform = view.camera.transform
            shading_filter.light_dir = transform.map(initial_light_dir)[:3]

    return wireframe_filter, shading_filter, attach_headlight


def make_states():
    shading = {'shading':None}, {'shading':'flat'}, {'shading':'smooth'}
    w, f = 'wireframe_only', 'faces_only'
    wireframe = {w:False, f:False}, {w:True, f:False}, {w:False, f:True}
    return shading, wireframe


def cycle_state(states, index):
    nindex = (index + 1) % len(states)
    return states[nindex], nindex


def add_shading_wireframe_interface(canvas, mesh, shading, wireframe):
    shading_states, wireframe_states = make_states()

    global shading_index
    global wireframe_index
    shading_index = shading_states.index({'shading':shading.shading})
    wireframe_index = wireframe_states.index({'wireframe_only':wireframe.wireframe_only, 'faces_only':wireframe.faces_only})

    @canvas.events.key_press.connect
    def on_key_press(event):
        global shading_index
        global wireframe_index
        if event.key == 's':
            state, shading_index = cycle_state(shading_states, shading_index)
            for attr, value in state.items():
                setattr(shading, attr, value)
            mesh.update()
        elif event.key == 'w':
            wireframe.enabled = not wireframe.enabled
            mesh.update()
        elif event.key == 'f':
            state, wireframe_index = cycle_state(wireframe_states, wireframe_index)
            for attr, value in state.items():
                setattr(wireframe, attr, value)
            mesh.update()
        elif event.key == '+' or event.key == '=':
            shading.shininess = min(100, shading.shininess + 10)
            mesh.update()
        elif event.key == '-':
            shading.shininess = max(0, shading.shininess - 10)
            mesh.update()


def is_interactive():
    """Identify if the current running Python is executing a script or running interactively

    Running a command from the command prompt, e.g., `python -c 'import brille.vis as bv; bv.something()`
    will be identified as an interactive session via this method.
    """
    import __main__ as main
    return not hasattr(main, '__file__')

