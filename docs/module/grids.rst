Grids
-----
There are likely myriad ways one could divide an :math:`N` dimensional
space but the two most basic building blocks are the parallelepiped
(or its :math:`N` dimensional equivalent) and the :math:`N`-simplex.
An :math:`N`-simplex always requires :math:`N+1` points
where a parallelepiped requires :math:`2^N`.

The parallelepiped lends itself to a regular grid while arrangements of `N`-simplexs can be regular or irregular.
In two dimensions the parallelepiped grid is a regular Cartesian grid:

.. tikz::
    :libs: calc

    \begin{tikzpicture}[
    dot/.style = {radius=0.04}, dotfill/.style={color=black, fill},
    ]
    \draw[step=1 cm] (-0.6,-0.2) grid (2.1,1.6);
    \coordinate (A) at (0,0);
    \coordinate (B) at (0,1);
    \coordinate (C) at (1,1);
    \coordinate (D) at (1,0);
    \coordinate (E) at (2,0);
    \coordinate (F) at (2,1);
    %
    \draw (A) -- (B) -- (C) -- (D) -- cycle;
    %
    \draw[dotfill] (A) circle [dot];
    \draw[dotfill] (B) circle [dot];
    \draw[dotfill] (C) circle [dot];
    \draw[dotfill] (D) circle [dot];
    \draw[dotfill] (E) circle [dot];
    \draw[dotfill] (F) circle [dot];
    %
    \end{tikzpicture}

while an irregular mesh of triangles is also possible:

.. tikz::
    :libs: calc

    \begin{tikzpicture}[%
    dot/.style = {radius=0.04},
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
    %
    \draw (A) -- (B) -- (C) -- (E) -- (F) -- (G) -- cycle;
    \draw (A) -- (D) -- (B);
    \draw (D) -- (C) -- (F) -- cycle;
    \draw (G) -- (D);
    \draw (B) -- (H) -- (C);
    \draw (H) -- (E);
    %
    \draw[dotfill] (A) circle [dot];
    \draw[dotfill] (B) circle [dot];
    \draw[dotfill] (C) circle [dot];
    \draw[dotfill] (D) circle [dot];
    \draw[dotfill] (E) circle [dot];
    \draw[dotfill] (F) circle [dot];
    \draw[dotfill] (G) circle [dot];
    \draw[dotfill] (H) circle [dot];
    %
    \end{tikzpicture}


The :py:mod:`brille._brille` module implements multiple grid types for linear interpolation
within the first or irreducible Brillouin zone.
In most cases the :py:class:`~brille._brille.BZTrellisQdc` should be used.

Since the grids are intended to be used with eigenvectors and their associated eigenvalues
each must support mixed real and complex data.
To allow this, each type of grid is exposed to Python three times with a suffix
`dd` for `(double, double)`, `dc` for `(double, std::complex<double>)`
and `cc` for `(std::complex<double>, std::complex<double>)`.

Regular grids
=============

Interpolation using grids constructed of regular parallelepipeds is possible,
but will lead to artefacts near the Brillouin zone boundary if the primitive
lattice has any non-orthogonal basis vectors.
This is such a large restriction that regular grids are no longer supported in :py:mod:`brille._brille`

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


Triangulated grids
==================

An alternative to the regular parallelepiped grid which may still have some utility is the
triangulation of a bounded space.
If the triangulated space is bounded by a polyhedron then the triangulation can always
represent it exactly, though possibly with sub-optimal triangulated cells.


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


There are two exposed implementations of the triangulated grid,
both detailed below under the headings `Simple`_ and `Hierarchy`_.

Simple
^^^^^^
The :py:class:`~brille._brille.BZMeshQdc` and its type-siblings implement a simple triangulated grid.
Locating the tetrahedron within the grid which contains a test point could require as many in-tetrahedron
checks as there are tetrahedra in the grid.
This class should be fine for use in applications where intra-grid-point interpolation is not required,
such as Brillouin zone integrations, but should be avoided when interpolation at random points is required.

.. autoclass:: brille._brille.BZMeshQdd

.. autoclass:: brille._brille.BZMeshQcc

.. autoclass:: brille._brille.BZMeshQdc
  :members:


Hierarchy
^^^^^^^^^
The :py:class:`~brille._brille.BZNestQdc` and its type-siblings implements a triangulated grid with
multiple triangulations of the Brillouin Zone of increasingly finer maximum tetrahedron size;
overlapping regions of these triangulations are identified and stored to simplify later finding
the smallest tetrahedron containing a given point.

Starting from the coarsest 'top' `layer` of tetrahedra, a containing tetrahedron should be quick to locate.
Then the list of overlapping next-finer `layer` tetrahedra can be used to locate the next containing tetrahedron.
This process repeats until the finest `layer` tetrahedron containing the point is found.
The total number of in-tetrahedron checks is then on the order of the sum of the average number of connected
tetrahedra at each layer over the number of layers; which should be chosen to be smaller than the total
number of tetrahedra at the finest layer.

.. Note::
   As implemented, the tetrahedra of two subsequent triangulations could have any overlapping relationship
   including one finer-tetrahedron intersecting with two (or more) coarser-tetrahedra.
   In practice it seems that the `TetGen` meshing algorithm always produces finer-tetrahedra which subdivide
   a coarser-tetrahedra, so that a finite-tree is formed between the `layers` of the `nest`.
   A class requiring this relationship should be able to extract memory and speed efficiencies which the
   current implementation sacrificed for greater relationsional flexibility.

If the :math:`i^\text{th}` layer has tetrahedra which on average contain :math:`\left<m\right>_i` overlapping 
number of tetrahedra in the :math:`(i+1)^\text{th}` layer, then the typical number of in-tetrahedra
checks required to find the finest-tetrahedra containing any point is

.. math::
   N_\text{checks} \propto \sum_{i=1}^n \left< m \right>_i

and the finest `layer` contains

.. math::
   N_\text{tetrahedra}^n = \prod_{i=1}^n \left< m \right>_i

tetrahedra. As long as :math:`N_\text{checks} < N_\text{tetrahedra}^n` this location method will be more efficient.

.. autoclass:: brille._brille.BZNestQdd

.. autoclass:: brille._brille.BZNestQcc

.. autoclass:: brille._brille.BZNestQdc
  :members:


.. _hybrid-grids:

Hybrid grids
============

A hybrid grid employs both a regular grid of parallelepipeds and, where the regular grid
passes the Brillouin zone boundary, individually triangulated cells.
This enables fast location of the cell containing any interpolation point and, if it
is triangulated, subsequent fast location of the containing tetrahedron.

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


The three exposed hybrid grid implementations have the same set of methods and are:

.. autoclass:: brille._brille.BZTrellisQdd

.. autoclass:: brille._brille.BZTrellisQcc

.. autoclass:: brille._brille.BZTrellisQdc
  :members:


Helper classes
==============

.. autoclass:: brille._brille.RotatesLike
