============
Introduction
============

Background
##########

Lattices
========

An infinite set of points which are indistinguishable for a specific set of
rigid translations form a lattice.

The easiest conceived lattice is one formed by the combinations of the whole
numbers on a set of Cartesian axes. In two dimensions, this is equivalent to
the intersections of lines on a piece of graphing paper.
If the vector pointing from any one such intersection point to its neighbour
along the :math:`x` axis is :math:`\mathbf{a}=[1\,0]`, and along the :math:`y`
axis is :math:`\mathbf{b}=[0\,1]`, then the allowed translations which result
in an indistinguishable lattice are
:math:`\mathbf{r} = i\mathbf{a} + j\mathbf{b}`
where both :math:`i` and :math:`j` are integers.

Although only translational invariance is necessary for a lattice, this
specific lattice has a number of additional invariant operations in the form of
rotations and rotoinversions that form a group,

.. math::

    \def \mat#1#2#3#4{\begin{pmatrix} #1 & #2\\#3 & #4\end{pmatrix}}
    \mathbb{G} = \mat{1}{0}{0}{1}, \mat{0}{-1}{1}{0}, \mat{0}{1}{1}{0}, \ldots

where for every element, :math:`G\in\mathbb{G}`, there exists an equivalent
product of elements of :math:`\mathbb{G}`. This allows a reduced representation
of :math:`\mathbb{G}` by a minimal set of its elements, called its generators,
which can be combined to produce all other elements.

Lattices like the one described above, where the repeated unit is point-like,
are called *Bravais* lattices and, in three dimensions, there are :math:`14`
unique Bravais lattices.

If the repeated unit is considered independent of the lattice, it is possible
to conceive of shapes which have their own rotation and rotoinversion
symmetries. There are :math:`32` unique sets of rotational symmetries in three
dimensions each of which forms a *Point group*.

Combining the symmetries of a Point group and a Bravais lattice gives rise to
a set of symmetry operations with is comprised of rotations, rotoinversions,
translations, screw axes, and glide planes which form a *Space group*, of
which there are :math:`230` distinct combinations in three dimensions.

All crystalline materials have a structure characterised by one of the
Space groups.


Dual lattices
-------------------

While the specific example of a :py:class:`~brille._brille.Lattice` given above
was described in terms of a physical lattice, Lattices are not limited to any
one space. For every Lattice there exists a dual lattice which is, effectively, its
inverse. While a direct lattice describes atom positions, its dual lattice
describes the relative orientations of *planes* of atoms.
Constructing such a dual lattice is straightforward following references easily
found elsewhere, but importantly where a lattice can be described by the lengths
of its basis vectors, the basis vectors of its dual lattice have units of
inverse, or reciprocal, length.
The :py:class:`~brille._brille.Lattice` class combines both the real space direct
lattice and its reciprocal space dual lattice into a single object.



Physical Properties
-------------------

The physical properties of any crystalline material have the same symmetry
as its Space group.
The goal of :py:mod:`brille` is to simplify leveraging the Space group symmetry
information in other software projects dealing with the physical properties
of crystalline materials.


Brillouin zone
==============

The region of space around each lattice point which is closer to that point than
to any other is its Wigner-Seitz cell.
As we can see for the case of a Lattice with hexagonal symmetry,
the Wigner-Seitz cells of all lattice points tile the space of the Lattice.

.. tikz::
    :libs: calc

    \begin{tikzpicture}[scale=5,dot/.style = {fill,radius=0.02},
    latpt/.style = {color=gray!50!white, fill},]
    \coordinate (astar) at (1,0);
    \coordinate (bstar) at (0.5,0.8660254037844386);
    \clip ($-0.1*(bstar)$) rectangle ($3*(astar)+2.1*(bstar)$);
    %
    \coordinate (C1) at ($0.3333*(astar)+0.3333*(bstar)$);
    \coordinate (C2) at ($-0.3333*(astar)+0.6667*(bstar)$);
    \coordinate (C3) at ($-0.6667*(astar)+0.3333*(bstar)$);
    \coordinate (C4) at ($-0.3333*(astar)-0.3333*(bstar)$);
    \coordinate (C5) at ($0.3333*(astar)-0.6667*(bstar)$);
    \coordinate (C6) at ($0.6667*(astar)-0.3333*(bstar)$);
    %
    \foreach \h in {-1,...,6}{
    \foreach \k in {-1,...,3}{
      \draw[latpt] ($\h*(astar) + \k*(bstar)$) circle[dot];
      \draw[color=yellow!75!black,dotted] ($\h*(astar)+\k*(bstar)$) +(C1)
        -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
    }}
    %
    \end{tikzpicture}

If we consider that the Lattice is a Reciprocal lattice, with basis vectors
:math:`\mathbf{a}^*` and :math:`\mathbf{b}^*`, then we can construct the
Wigner-Seitz cell for any lattice point
:math:`\boldsymbol{\tau} = h\mathbf{a}^* + k\mathbf{b}^*`
by finding the intersections of all planes, each described by
the point :math:`\boldsymbol{\tau}+\frac{i}{2}\mathbf{a}^*+\frac{j}{2}\mathbf{b}^*`
and normal vector :math:`i\mathbf{a}^*+j\mathbf{b}^*` with :math:`i,j` integers;
the Wigner-Seitz cell is bounded by all such planes *without* any planes closer
to :math:`\boldsymbol{\tau}`.

This is also the definition of the first Brillouin zone, where higher-order
Brillouin zones are the regions between planes successively further from
:math:`\mathbf{G}`.

.. tikz::
    :libs: calc

    \begin{tikzpicture}[scale=5,dot/.style = {fill,radius=0.02},
    latpt/.style = {color=gray!50!white, fill},]
    \coordinate (astar) at (1,0);
    \coordinate (bstar) at (0.5,0.8660254037844386);
    \clip ($-0.1*(bstar)$) rectangle ($3*(astar)+2.1*(bstar)$);
    %
    \coordinate (C1) at ($0.3333*(astar)+0.3333*(bstar)$);
    \coordinate (C2) at ($-0.3333*(astar)+0.6667*(bstar)$);
    \coordinate (C3) at ($-0.6667*(astar)+0.3333*(bstar)$);
    \coordinate (C4) at ($-0.3333*(astar)-0.3333*(bstar)$);
    \coordinate (C5) at ($0.3333*(astar)-0.6667*(bstar)$);
    \coordinate (C6) at ($0.6667*(astar)-0.3333*(bstar)$);
    %
    \foreach \h in {-1,...,6}{
    \foreach \k in {-1,...,3}{
      \draw[latpt] ($\h*(astar) + \k*(bstar)$) circle[dot];
      \draw[color=yellow!75!black,dotted] ($\h*(astar)+\k*(bstar)$) +(C1)
        -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
    }}
    %
    \draw[<->,very thick, gray] (bstar) -- (0,0) node[near start,above left]{$b^*$}  -- (astar) node[near end,above]{$a^*$};
    %
    \coordinate (G) at ($(astar)+(bstar)$);
    \draw[->, very thick, dashed] (0,0) -- (G) node[midway,above left] {$\tau$};
    \draw[color=yellow!75!black, line width=1mm] (G) +(C1) -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
    \end{tikzpicture}

Since the properties of the lattice follow the periodicity of the lattice, any
measurable quantity must repeat from one first Brillouin zone to the next.
This allows for descriptions of the physical properties which depend on, e.g.,
a reduced momentum transfer :math:`\mathbf{q} = \mathbf{Q}-\mathbf{G}`

.. tikz::
    :libs: calc

    \begin{tikzpicture}[scale=5,dot/.style = {fill,radius=0.02},
    latpt/.style = {color=gray!50!white, fill},]
    \coordinate (astar) at (1,0);
    \coordinate (bstar) at (0.5,0.8660254037844386);
    \clip ($-0.1*(bstar)$) rectangle ($3*(astar)+2.1*(bstar)$);
    %
    \coordinate (C1) at ($0.3333*(astar)+0.3333*(bstar)$);
    \coordinate (C2) at ($-0.3333*(astar)+0.6667*(bstar)$);
    \coordinate (C3) at ($-0.6667*(astar)+0.3333*(bstar)$);
    \coordinate (C4) at ($-0.3333*(astar)-0.3333*(bstar)$);
    \coordinate (C5) at ($0.3333*(astar)-0.6667*(bstar)$);
    \coordinate (C6) at ($0.6667*(astar)-0.3333*(bstar)$);
    %
    \foreach \h in {-1,...,6}{
    \foreach \k in {-1,...,3}{
      \draw[latpt] ($\h*(astar) + \k*(bstar)$) circle[dot];
      \draw[color=yellow!75!black,dotted] ($\h*(astar)+\k*(bstar)$) +(C1)
        -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
    }}
    %
    \draw[<->,very thick, gray] (bstar) -- (0,0) node[near start,above left]{$b^*$}  -- (astar) node[near end,above]{$a^*$};
    %
    \coordinate (G) at ($(astar)+(bstar)$);
    \draw[->, very thick, dashed] (0,0) -- (G) node[midway,above left] {$\tau$};
    \draw[color=yellow!75!black, line width=1mm] (G) +(C1) -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
    \coordinate (q) at ($0.1*(astar)-0.3*(bstar)$);
    \coordinate (Q) at ($(G)+(q)$);
    \draw[->,very thick,color=black] (0,0) -- (Q) node[midway,below right] {$Q$};
    \draw[->,very thick,color=black] (G) -- (Q) node[midway,below right] {$q$};
    \end{tikzpicture}


Irreducible first Brillouin zone
--------------------------------

The first Brillouin zone may contain redundant information depending on
the Point group symmetry of the Space group.
If the two-dimensional hexagonal lattice above possesses a six-fold rotation
axis perpendicular to the plane, so that the information within each first
Brillouin zone is repeated six times, then the zone can be *reduced*.
The definition of an irreducible zone is not unique, but one choice for
this Reciprocal lattice is shown below

.. tikz::
    :libs: calc

    \begin{tikzpicture}[scale=5,dot/.style = {fill,radius=0.02},
    latpt/.style = {color=gray!50!white, fill},]
    \coordinate (astar) at (1,0);
    \coordinate (bstar) at (0.5,0.8660254037844386);
    \clip ($-0.1*(bstar)$) rectangle ($3*(astar)+2.1*(bstar)$);
    %
    \coordinate (C1) at ($0.3333*(astar)+0.3333*(bstar)$);
    \coordinate (C2) at ($-0.3333*(astar)+0.6667*(bstar)$);
    \coordinate (C3) at ($-0.6667*(astar)+0.3333*(bstar)$);
    \coordinate (C4) at ($-0.3333*(astar)-0.3333*(bstar)$);
    \coordinate (C5) at ($0.3333*(astar)-0.6667*(bstar)$);
    \coordinate (C6) at ($0.6667*(astar)-0.3333*(bstar)$);
    %
    \foreach \h in {-1,...,6}{
    \foreach \k in {-1,...,3}{
      \draw[latpt] ($\h*(astar) + \k*(bstar)$) circle[dot];
      \draw[color=yellow!75!black,dotted] ($\h*(astar)+\k*(bstar)$) +(C1)
        -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
    }}
    %
    \draw[<->,very thick, gray] (bstar) -- (0,0) node[near start,above left]{$b^*$}  -- (astar) node[near end,above]{$a^*$};
    %
    \coordinate (G) at ($(astar)+(bstar)$);
    \draw[->, very thick, dashed] (0,0) -- (G) node[midway,above left] {$\tau$};
    \draw[color=yellow!75!black,line width=1mm] (G) -- +(C4) -- +(C5) -- cycle;
    %\draw[color=yellow!75!black,dashed]         (G) -- +(C1) -- +(C2) -- cycle;
    %\draw[color=yellow!75!black,dashed]         (G) -- +(C2) -- +(C3) -- cycle;
    %\draw[color=yellow!75!black,dashed]         (G) -- +(C3) -- +(C4) -- cycle;
    %\draw[color=yellow!75!black,dashed]         (G) -- +(C5) -- +(C6) -- cycle;
    %\draw[color=yellow!75!black,dashed]         (G) -- +(C6) -- +(C1) -- cycle;
    \coordinate (q) at ($0.1*(astar)-0.3*(bstar)$);
    \coordinate (Q) at ($(G)+(q)$);
    \draw[->,very thick,color=black] (0,0) -- (Q) node[midway,below right] {$Q$};
    \draw[->,very thick,color=black] (G) -- (Q) node[midway,below right] {$q$};
    \end{tikzpicture}


The properties at an arbitrary momentum transfer :math:`\mathbf{Q}` can be
related to those within the irreducible first Brillouin zone by

.. math::

    \mathbf{Q} = G \mathbf{q}_\text{ir} + \boldsymbol{\tau}

where :math:`G` is one of the Point group operators,
:math:`\mathbf{q}_\text{ir}` is a vector within the irreducible first
Brillouin zone, and :math:`\boldsymbol{\tau}` is a Reciprocal lattice point.

Since the irreducible first Brillouin zone contains all of the information about
the physical properties of a material, it can and should be used by projects
aiming to model those properties efficiently.
To help in this task, :py:mod:`brille` defines :py:class:`~brille._brille.BrillouinZone` to
construct the first Brillouin zone and an irreducible Brillouin zone for any
Reciprocal lattice.


Inelastic Neutron Scattering
----------------------------

Inelastic neutron scattering is an experimental technique which measures the
probability of transitions between states of a condensed matter system, which
in turn can tell us about the types and strengths of interactions within the
material.

Inelastic neutron scattering benefits from the use of neutrons with wavelengths
comparable to typical interatomic spacings *and* energies comparable to typical
energy levels of condensed matter systems.

The straightforward comparison of intensity measured on a neutron spectrometer
and favourable wavelength and energy of available neutrons compensates for the
difficulty of neutron production compared to, e.g., x-rays which are easier to
produce but can not have both favourable wavelengths and energies in the same
photon.

The difficulty of producing neutron beams led to the development of instruments
like the Direct Geometry Time of Flight neutron spectrometer. Such instruments
have an array of detectors at fixed positions and detect changes in the
neutron energy by measuring the time it takes for a detected neutron to arrive
at the detector. By knowing the neutron's initial, :math:`\mathbf{k}_\text{i}`,
and final momentum, :math:`\mathbf{k}_\text{f}` it is
straightforward to work out the momentum and energy transferred to the sample.

.. math::

    \begin{aligned}
    \mathbf{Q} & = \mathbf{k}_\text{i} - \mathbf{k}_\text{f} \\
    E & = \frac{\hbar^2}{2m_\text{n}}\left(k_\text{i}^2 - k_\text{f}^2\right)
    \end{aligned}


Motivation
##########

Through the use of one or more choppers, Direct Geometry Time of Flight
spectrometers select a single :math:`\mathbf{k}_\text{i}` for all neutrons which
interact with the sample before being counted in a detector. Each detector is
at a unique set of spherical angles :math:`(\theta,\phi)` relative to
:math:`\hat{\mathbf{k}_\text{i}}` and therefore each counts neutrons with
a unique :math:`\hat{\mathbf{k}}_\text{f}`. As a result each detector measures
along a path through reciprocal :math:`(\mathbf{Q},E)` space which is
constrained by the kinematic relations listed above.

Theoretical models of interactions in a material typically involve solving
an eigenvalue problem for a given :math:`\mathbf{Q}` and are therefore best
suited for simulating along :math:`(\mathbf{Q},E)` paths with
constant-:math:`\mathbf{Q}`.
:py:mod:`brille` aims to help such models by reducing the number of
:math:`\mathbf{Q}` points where they must perform their (typically expensive)
calculation and interpolates their results onto the :math:`(\mathbf{Q},E)` paths
measured during experiments.
To accomplish this, a number of polyhedron-filling connected grids are defined;
notably :py:class:`~brille._brille.BZTrellisQdc` and similar variants.
The model calculation is evaluated for :math:`\mathbf{Q}` points defined by
the vertices of the polyhedra which comprise the grid.
Interpolation is done then by finding the polyhedron which encloses the desired
point and linearly weighting the pre-calculated values at the vertices by the
distance from the desired point to each vertex.


Irreducible Brillouin zone interpolation
========================================

Since physical properties of crystalline solids are unique only within the 
first irreducible Brillouin zone, only :math:`\mathbf{Q}` points within this
region need be calculated by the expensive model calculation. Thus, the 
polyhedron-filling connected grids takes a first Brillouin zone
polyhedron or an irreducible Brillouin zone polyhedron plus, e.g., a maximum
distance between grid nodes or a maximum grid cell volume, and define a grid.

Cartesian grid
--------------

.. tikz::
    :libs: calc

    \begin{tikzpicture}[
    dot/.style = {radius=0.04}, dotfill/.style={color=black, fill},
    adot/.style = {color=red},
    bdot/.style = {color=blue},
    cdot/.style = {color=green!50!black},
    ddot/.style = {color=orange},
    a/.style = {fill=red!50!white},
    b/.style = {fill=blue!50!white},
    c/.style = {fill=green!50!black!50!white},
    d/.style = {fill=orange!50!white},
    ]
    \draw[step=1 cm] (-0.6,-0.2) grid (2.1,1.6);
    \coordinate (A) at (0,0);
    \coordinate (B) at (0,1);
    \coordinate (C) at (1,1);
    \coordinate (D) at (1,0);
    \coordinate (E) at (2,0);
    \coordinate (F) at (2,1);
    \coordinate (X) at (0.6, 0);
    \coordinate (Y) at (0, 0.2);
    \coordinate (XY) at ($(X)+(Y)$);
    %
    \fill[a] (C) rectangle (XY);
    \fill[b] (D) rectangle (XY);
    \fill[c] (A) rectangle (XY);
    \fill[d] (B) rectangle (XY);
    %
    \draw (A) -- (B) -- (C) -- (D) -- cycle;
    %
    \draw[adot,a] (A) circle [dot];
    \draw[bdot,b] (B) circle [dot];
    \draw[cdot,c] (C) circle [dot];
    \draw[ddot,d] (D) circle [dot];
    \draw[color=black,fill=white] (XY) circle [dot];
    \draw[color=black, fill] (E) circle [dot];
    \draw[color=black, fill] (F) circle [dot];
    %
    \end{tikzpicture}

One simple approach to defining a grid within a polyhedron is to

#. define one vertex of the polyhedron as the origin,
#. find the vertex farthest away from the origin
#. subdivide the rectangular prism defined by these two points.

Such a grid has the advantage that for all space within it, the closest grid
point(s) can be calculated analytically. This lends itself to fast neighbour
location and fast linear interpolation.

A disadvantage to such a grid is that it can only be commensurate with
polyhedra which are also rectangular prisms, which is the case only for
irreducible first Brillouin zones of primitive cubic, primitive tetragonal, 
and primitive orthorhombic space groups.
When the grid is not commensurate with the polyhedron it is likely to introduce
unmanageable artifacts in any interpolation result.

The disadvantages of the basic Cartesian grid are so restrictive that
:py:mod:`brille` does not implement a three-dimensional Cartesian grid object.

:math:`n`-simplex grid
----------------------

.. tikz::
    :libs: calc

    \begin{tikzpicture}[%
    dot/.style = {radius=0.04},
    afill/.style={color=purple!50!white},
    bfill/.style={color=yellow!50!green!50!white},
    cfill/.style={color=teal!50!white},
    adot/.style={color=purple,fill=purple!50!white},
    bdot/.style={color=yellow!50!green,fill=yellow!50!green!50!white},
    cdot/.style={color=teal,fill=teal!50!white},
    dotfill/.style={color=black, fill},
    ]
    \coordinate (A) at (0,0);
    \coordinate (B) at (1,-0.2);
    \coordinate (C) at (2,1);
    \coordinate (D) at (0.7,0.9);
    \coordinate (E) at ($(C)+(0.5,0.8)$);
    \coordinate (F) at ($(D)+(0.3,0.8)$);
    \coordinate (G) at ($(A)+(-0.3,1.2)$);
    \coordinate (H) at ($(B)+(1.4,0.3)$);
    \coordinate (XY) at ($0.3*(B)+0.2*(C)+0.5*(D)$);
    %
    \fill[afill] (C) -- (D) -- (XY) -- cycle;
    \fill[bfill] (B) -- (XY) -- (D) -- cycle;
    \fill[cfill] (B) -- (C) -- (XY) -- cycle;
    %
    \draw (A) -- (B) -- (C) -- (E) -- (F) -- (G) -- cycle;
    \draw (A) -- (D) -- (B);
    \draw (D) -- (C) -- (F) -- cycle;
    \draw (G) -- (D);
    \draw (B) -- (H) -- (C);
    \draw (H) -- (E);
    %
    \draw[dotfill] (A) circle [dot];
    \draw[adot] (B) circle [dot];
    \draw[bdot] (C) circle [dot];
    \draw[cdot] (D) circle [dot];
    \draw[color=black,fill=white] (XY) circle [dot];
    \draw[dotfill] (E) circle [dot];
    \draw[dotfill] (F) circle [dot];
    \draw[dotfill] (G) circle [dot];
    \draw[dotfill] (H) circle [dot];
    %
    \end{tikzpicture}

Another straightforward approach to defining a grid within a polyhedron is the
use of a tetrahedral tiling. Tetrahedra being the three dimensional simplex.
Creating such a tiling with nice properties is nontrivial, so :py:mod:`brille` uses the
`TetGen <http://tetgen.org>`_ library to do the heavy lifting.

Tetrahedral tilings have the advantage that they can be made commensurate with
any polyhedron, and therefore never introduce unmanageable artefacts when
interpolating near their surfaces.
But they lack the ability to calculate which tetrahedron contains a specified
point.
So interpolating with a tetrahedral tiling is either slow or requires
substantial meta-information to be determined in advance.

The classes :py:class:`~brille._brille.BZMeshQdc` and :py:class:`~brille._brille.BZNestQdc` 
implement 3-D :math:`n`-simplex grids which fills and does not extend beyond
the boundaries of a [irreducible] Brillouin zone. The differences between the
two classes relate to how their (meta)data is stored, either in a flat or tree
format as described in the :doc:`grids page <module/grids>`.


Hybrid grid
-----------

.. tikz::
    :libs: calc

    \begin{tikzpicture}[scale=3,%
    	dot/.style = {radius=0.0133},
     	cell/.style = {color=black},
      s0/.style = {color=red, fill=red!50!white},
      s1/.style = {color=blue, fill=blue!50!white},
      s2/.style = {color=green!50!black, fill=green!50!black!50!white},
      s3/.style = {color=orange, fill=orange!50!white},
      c0/.style = {color=purple, fill=purple!50!white},
      c1/.style = {color=yellow!50!green, fill=yellow!50!green!50!white},
      c2/.style = {color=teal, fill=teal!50!white},
    ]
    \draw[step=2.8867513459mm, color=black!30!white, dashed, very thin] (0,0) grid (0.9,1.05);
    \foreach \i in {0,...,3} {\foreach \j in {0,...,4} {\coordinate (g\i\j) at (0.28867513459*\i, 0.28867513459*\j);}}
    % BZ boundary points
    \coordinate (BZ1) at (30:1);
    \coordinate (BZ0) at ($(BZ1) +(0,-0.5)$);
    \coordinate (BZ2) at (90:1);
    % Brillouin zone boundary
    \draw[color=yellow!75!black, line width=1mm] (BZ0) -- (BZ1) -- (BZ2);
    % regular cell points
    \coordinate (dx) at (2.8867513459mm, 0);
    \coordinate (dy) at (0, 2.8867513459mm);
    % extra triangulation points
    \coordinate (e0) at (BZ1);
    \coordinate (e1) at ($0.54*(g32)+0.46*(g22)$);
    \coordinate (e2) at ($(BZ1)+(150:0.333333)$);
    \coordinate (e3) at ($(BZ1)+(150:0.666667)$);
    \coordinate (e4) at ($0.80*(g13)+0.20*(g03)$);
    \coordinate (e5) at (BZ2);
    % 'regular' interpolation points
    \coordinate (r0) at (g10);
    \coordinate (r1) at (g20);
    \coordinate (r2) at (g21);
    \coordinate (r3) at (g11);
    % 'simplex' interpolation points
    \coordinate (t0) at (g12);
    \coordinate (t1) at (e2);
    \coordinate (t2) at (e3);
    % x1
    \coordinate (x1) at ($0.25*(t0)+0.33*(t1)+0.42*(t2)$);
    \fill[c2] (t0) -- (t1) -- (x1) -- cycle;
    \fill[c0] (t1) -- (t2) -- (x1) -- cycle;
    \fill[c1] (t2) -- (t0) -- (x1) -- cycle;
    % x2
    \coordinate (x2) at ($0.15*(r0) + 0.39*(r1) + 0.06*(r2) + 0.4*(r3)$);
    \fill[s0] (r2) rectangle (x2);
    \fill[s1] (r3) rectangle (x2);
    \fill[s2] (r0) rectangle (x2);
    \fill[s3] (r1) rectangle (x2);
    % full cells
    \foreach \pt in {(g00), (g10), (g20), (g01), (g11)} { \draw[cell] \pt rectangle +(g11); }
    % triangulated cells
    \draw [cell] (e1) -- (e0) -- (g31) -- (e1) -- (g21);
    \draw [cell] (e1) -- (g22) -- (e2) -- cycle;
    \draw [cell] (e2) -- (g22) -- (g12) -- (e2) -- (e3) -- (g12) -- cycle;
    \draw [cell] (g03) -- (g02) -- (e3) -- (g03) -- (e4) -- (e3);
    \draw [cell] (e4) -- (g03) -- (e5) -- cycle;
    % mesh/grid points
    \foreach \pt in {(x1), (x2)} {\draw [fill=white] \pt circle[dot];}
    \foreach \i in {0,...,2} {\draw[c\i] (t\i) circle[dot];}
    \foreach \i in {0,...,3} {\draw[s\i] (r\i) circle[dot];}
    \foreach \pt in {(g00), (g30), (g31), (g01), (g02), (g22), (g03), (e0), (e1), (e4), (e5)} {\draw[fill=lightgray] \pt circle[dot];}
    \end{tikzpicture}

An alternative approach is to combine a Cartesian grid with a :math:`n`-simplex
grid. Such a grid has its rectangular-prism cells replaced by triangulated
truncated-rectangular-prisms on the surface of the polyhedron.

Such a construction has the advantage of direct calculation of the cell which
contains any given point with a much-faster search over only those tetrahedra
within the cell if the rectangular-prism passes the surface of the polyhedron.

The class :py:class:`~brille._brille.BZTrellisQdc` implements a hybrid grid in three dimensions
which fills and does not extend beyond the boundaries of a [irreducible]
Brillouin zone.
