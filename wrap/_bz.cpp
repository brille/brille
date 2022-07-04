/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */
#include <pybind11/pybind11.h>
#include "bz.hpp"
#include "lattice_dual.hpp"
#include "_c_to_python.hpp"
#include "_array.hpp"
#include "array_l_.hpp"

namespace py = pybind11;

void wrap_brillouinzone(py::module & m){
  using namespace pybind11::literals; // bring in "[name]"_a to be interpreted as py::arg("[name]")
  using namespace brille;
  using namespace brille::lattice;
  using CLS = BrillouinZone;
  py::class_<CLS> cls(m,"BrillouinZone",R"pbdoc(
    Construct and hold a first Brillouin zone and, optionally and by default,
    an irreducible Brillouin zone.

    The region closer to a given lattice point than to any other is the
    Wigner-Seitz cell of that lattice. The same construction is one possible
    first Brillouin zone of a reciprocal lattice and is used within ``brille``.
    For example, a two-dimensional hexagonal lattice has a first Brillouin
    zone which is a hexagon:

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
        \coordinate (G) at ($(astar)+(bstar)$);
        \draw[color=yellow!75!black, line width=1mm] (G) +(C1) -- +(C2) -- +(C3) -- +(C4) -- +(C5) -- +(C6) -- cycle;
        \end{tikzpicture}


    Since all physical properties of a crystal must have the same periodicity
    as its lattice, the powerful feature of the first Brillouin zone is that it
    encompasses a region of reciprocal space which must fully represent all
    of reciprocal space.

    Most crystals contain rotational or rotoinversion symmetries in addition to
    the translational ones which give rise to the first Brillouin zone. These
    symmetries are the pointgroup of the lattice and enforce that the properties
    of the crystal also have the same symmetry. The first Brillouin zone,
    therefore, typically contains redundant information.

    An irreducible Brillouin zone is a subsection of the first Brillouin zone
    which contains the minimal part required to have only unique crystal
    properties. This class can find an irreducible Brillouin zone for any
    crystal lattice. In the example of the hexagonal lattice there are six
    equivalent irreducible Brillouin zones one of which is:

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
        \coordinate (G) at ($(astar)+(bstar)$);
        \draw[color=yellow!75!black,line width=1mm] (G) -- +(C4) -- +(C5) -- cycle;
        \end{tikzpicture}

    Parameters
    ----------
    lattice: :py:class:`brille._brille.Reciprocal`
        The reciprocal space lattice for which a Brillouin zone will be found
    use_primitive: bool
        If the provided :py:class:`brille._brille.Reciprocal` lattice is a
        conventional Bravais lattice, this parameter controls whether the
        equivalent primitive Bravais lattice should be used to find the first
        Brillouin zone. This is ``True`` by default and should only be modified
        for testing purposes.
    search_length: int
        The Wigner-Seitz construction of the first Brillouin zone finds the
        volume of space closer to a chosen reciprocal lattice point than any
        other reciprocal lattice point. This is accomplished by successively
        dividing the space by planes halfway between the chosen point and a
        subset of all other planes. The subset used is controlled by
        `search_length` and is every unique :math:`(\pm s_i\,0\,0)`,
        :math:`(0\,\pm s_j\,0)`, :math:`(0\,0\,\pm s_k)`,
        :math:`(\pm s_i\,\pm s_j\,0)`, :math:`(\pm s_i\,0\,\pm s_k)`,
        :math:`(0\,\pm s_j\,\pm s_k)`, :math:`(\pm s_i\,\pm s_j\,\pm s_k)` for
        :math:`1 \le s_\alpha \le` `search_length`.
        If the reciprocal lattice is primitive then the default `search_length`
        of ``1`` should always give the correct first Brillouin zone.
        For extra assurance that the correct first Brillouin zone is found, the
        procedure is internally repeated with `search_length` incremented by
        one and an error is raised if the two constructed polyhedra have
        different volumes.
    time_reversal_symmetry: bool
        Controls whether time reversal symmetry should be added to pointgroups
        lacking space inversion. This affects the found irreducible Brillouin
        zone for such systems. To avoid inadvertently adding time reversal
        symmetry when it is not appropriate, this is ``False`` by default.
    wedge_search: bool
        Controls whether an irreducible Brillouin zone should be found. With
        this set to ``False`` the returned :py:class:`brille._brille.BrillouinZone`
        will only contain the first Brillouin zone. If ``True`` the pointgroup
        symmetry operations will be used to identify *an* irreducible Brillouin
        zone as well. If the provided lattice's parameters do not match the
        symmetry of the pointgroup (e.g., a lattice which should be tetragonal
        like :math:`I4/mmm` but constructed with :math:`\gamma=120^\circ`) the
        algorithm will fail to find an appropriate irreducible Brillouin zone
        and an error will be raised. (Set to ``True`` by default).
  )pbdoc");
  cls.def(py::init([](
    const lattice::Lattice<double> & lat,
    const bool use_primitive,
    const int search_length,
    const bool time_reversal_symmetry,
    const bool wedge_search,
    const bool divide_primitive
    ){
    auto cfg = BrillouinZoneConfig();
    cfg.primitive(use_primitive);
    cfg.divide_extent(search_length);
    cfg.divide_primitive(divide_primitive);
    cfg.time_reversal(time_reversal_symmetry);
    cfg.wedge_search(wedge_search);
    return BrillouinZone(lat, cfg);
  }),
  "lattice"_a,
  "use_primitive"_a=true,
  "search_length"_a=1,
  "time_reversal_symmetry"_a=false,
  "wedge_search"_a=true,
  "divide_primitive"_a=true
  );
  cls.def(py::init([](
              const lattice::Lattice<double> & lat,
              const approx_float::Config ac,
              const bool use_primitive,
              const int search_length,
              const bool time_reversal_symmetry,
              const bool wedge_search,
              const bool divide_primitive
          ){
            auto cfg = BrillouinZoneConfig();
            cfg.primitive(use_primitive);
            cfg.divide_extent(search_length);
            cfg.divide_primitive(divide_primitive);
            cfg.time_reversal(time_reversal_symmetry);
            cfg.wedge_search(wedge_search);
            return BrillouinZone(lat, cfg, ac);
          }),
          "lattice"_a,
          "approx_config"_a,
          "use_primitive"_a=true,
          "search_length"_a=1,
          "time_reversal_symmetry"_a=false,
          "wedge_search"_a=true,
          "divide_primitive"_a=true
  );

  // return the internal Lattice object
  cls.def_property_readonly("lattice", [](const CLS &b){ return b.get_lattice();},
  R"pbdoc(
  Returns the defining :py:class:`brille._brille.Lattice` lattice
  )pbdoc");

  // access the polyhedra directly
  cls.def_property_readonly("polyhedron",&CLS::get_polyhedron,
  R"pbdoc(
  Returns the first Brillouin zone :py:class:`brille._brille.Polyhedron`
  )pbdoc");
  cls.def_property_readonly("ir_polyhedron",[](const CLS &b){return b.get_ir_polyhedron(true);},
  R"pbdoc(
  Returns the irreducible Brillouin zone :py:class:`brille._brille.Polyhedron`

  Returns
  -------
  :py:class:`brille._brille.Polyhedron`
    If no irreducible Brillouin zone was requested at construction, the returned
    polyhedron is that of the first Brillouin zone instead.
  )pbdoc");
  cls.def_property_readonly("ir_polyhedron_generated",[](const CLS &b){return b.get_ir_polyhedron(false);},
  R"pbdoc(
  Returns the found irreducible Brillouin zone :py:class:`brille._brille.Polyhedron`

  If the lattice pointgroup does not contain the space inversion operator
  the internally held 'irreducible' polyhedron is only half of the real
  irreducible polyhedron. This method gives access to the polyhedron found by
  the algorithm before being doubled for output.
  )pbdoc");

  // first Brillouin zone polyhedron
  cls.def_property_readonly("normals",
  [](const CLS &b){return brille::a2py(b.get_normals().hkl());},
  R"pbdoc(
  Return the first Brillouin zone face normals in rlu
  )pbdoc");
  cls.def_property_readonly("normals_invA",
  [](const CLS &b){return brille::a2py(b.get_normals().xyz());},
  R"pbdoc(
  Return the first Brillouin zone face normals in inverse ångstrom
  )pbdoc");
  cls.def_property_readonly("normals_primitive",
  [](const CLS &b){return brille::a2py(b.get_primitive_normals().hkl());},
  R"pbdoc(
  Return the first Brillouin zone face normals in primitive-lattice rlu
  )pbdoc");
  cls.def_property_readonly("points",
  [](const CLS &b){return brille::a2py(b.get_points().hkl());},
  R"pbdoc(
  Return the first Brillouin zone face centres in rlu
  )pbdoc");
  cls.def_property_readonly("points_invA",
  [](const CLS &b){return brille::a2py(b.get_points().xyz());},
  R"pbdoc(
  Return the first Brillouin zone face centres in inverse ångstrom
  )pbdoc");
  cls.def_property_readonly("points_primitive",
  [](const CLS &b){return brille::a2py(b.get_primitive_points().hkl());},
  R"pbdoc(
  Return the first Brillouin zone face centres in primitive-lattice rlu
  )pbdoc");
  cls.def_property_readonly("vertices",
  [](const CLS &b){return brille::a2py(b.get_vertices().hkl());},
  R"pbdoc(
  Return the first Brillouin zone unique face corners in rlu
  )pbdoc");
  cls.def_property_readonly("vertices_invA",
  [](const CLS &b){return brille::a2py(b.get_vertices().xyz());},
  R"pbdoc(
  Return the first Brillouin zone unique face corners in inverse ångstrom
  )pbdoc");
  cls.def_property_readonly("half_edge_points",
  [](const CLS &b){return brille::a2py(b.get_half_edges().hkl());},
  R"pbdoc(
  Return the first Brillouin zone face edge centres in rlu
  )pbdoc");
  cls.def_property_readonly("half_edge_points_invA",
  [](const CLS &b){return brille::a2py(b.get_half_edges().xyz());},
  R"pbdoc(
  Return the first Brillouin zone face edge centres in inverse ångstrom
  )pbdoc");
  cls.def_property_readonly("vertices_primitive",
  [](const CLS &b){return brille::a2py(b.get_primitive_vertices().hkl());},
  R"pbdoc(
  Return the first Brillouin zone unique face corners in primitive-lattice rlu
  )pbdoc");
  cls.def_property_readonly("faces_per_vertex",
  &CLS::get_faces_per_vertex,
  R"pbdoc(
  Return the first Brillouin zone face indices for each unique face corner
  )pbdoc");
  cls.def_property_readonly("vertices_per_face",
  &CLS::get_vertices_per_face,
  R"pbdoc(
  Return the first Brillouin zone face corner indices for each face
  )pbdoc");

  // irreducible first Brillouin zone polyhedron
  cls.def_property_readonly("ir_normals",
  [](const CLS &b){return brille::a2py(b.get_ir_normals().hkl());},
  R"pbdoc(
  Return the irreducible Brillouin zone face normals in rlu
  )pbdoc");
  cls.def_property_readonly("ir_normals_invA",
  [](const CLS &b){return brille::a2py(b.get_ir_normals().xyz());},
  R"pbdoc(
  Return the irreducible Brillouin zone face normals in inverse ångstrom
  )pbdoc");
  cls.def_property_readonly("ir_normals_primitive",
  [](const CLS &b){return brille::a2py(b.get_ir_primitive_normals().hkl());},
  R"pbdoc(
  Return the irreducible Brillouin zone face normals in primitive-lattice rlu
  )pbdoc");
  cls.def_property_readonly("ir_points",
  [](const CLS &b){return brille::a2py(b.get_ir_points().hkl());},
  R"pbdoc(
  Return the irreducible Brillouin zone face centres in rlu
  )pbdoc");
  cls.def_property_readonly("ir_points_invA",
  [](const CLS &b){return brille::a2py(b.get_ir_points().xyz());},
  R"pbdoc(
  Return the irreducible Brillouin zone face centres in inverse ångstrom
  )pbdoc");
  cls.def_property_readonly("ir_points_primitive",
  [](const CLS &b){return brille::a2py(b.get_ir_primitive_points().hkl());},
  R"pbdoc(
  Return the irreducible Brillouin zone face centres in primitive-lattice rlu
  )pbdoc");
  cls.def_property_readonly("ir_vertices",
  [](const CLS &b){return brille::a2py(b.get_ir_vertices().hkl());},
  R"pbdoc(
  Return the irreducible Brillouin zone unique face corners in rlu
  )pbdoc");
  cls.def_property_readonly("ir_vertices_invA",
  [](const CLS &b){return brille::a2py(b.get_ir_vertices().xyz());},
  R"pbdoc(
  Return the irreducible Brillouin zone unique face corners in inverse ångstrom
  )pbdoc");
  cls.def_property_readonly("ir_vertices_primitive",
  [](const CLS &b){return brille::a2py(b.get_ir_primitive_vertices().hkl());},
  R"pbdoc(
  Return the irreducible Brillouin zone unique face corners in primitive-lattice rlu
  )pbdoc");
  cls.def_property_readonly("ir_faces_per_vertex",
  &CLS::get_ir_faces_per_vertex,
  R"pbdoc(
  Return the irreducible Brillouin zone face index per unique face corner
  )pbdoc");
  cls.def_property_readonly("ir_vertices_per_face",
  &CLS::get_ir_vertices_per_face,
  R"pbdoc(
  Return the irreducible Brillouin zone unique face corners per face
  )pbdoc");

  // irreducible reciprocal space wedge
  cls.def_property_readonly("wedge_normals",
  [](const CLS &b){return brille::a2py(b.get_ir_wedge_normals().hkl());},
  R"pbdoc(
  Return the normals of the irreducible wedge rlu
  )pbdoc");
  cls.def_property_readonly("wedge_normals_invA",
  [](const CLS &b){return brille::a2py(b.get_ir_wedge_normals().xyz());},
  R"pbdoc(
  Return the normals of the irreducible wedge inverse ångstrom
  )pbdoc");
  cls.def_property_readonly("wedge_normals_primitive",
  [](const CLS &b){return brille::a2py(b.get_primitive_ir_wedge_normals().hkl());},
  R"pbdoc(
  Return the normals of the irreducible wedge primitive-lattice rlu
  )pbdoc");

  // check whether one or more points are inside
  cls.def("isinside",[](CLS &b, py::array_t<double> p){
    // brille::Array<double> sp = brille::py2a(p);
    brille::Array2<double> sp = brille::py2a2(p);
    if (sp.size(sp.ndim()-1) != 3)
      throw std::runtime_error("The last dimension must have size 3");
    auto pv = lattice::LVec<double>(LengthUnit::inverse_angstrom, b.get_lattice(), sp); // no copy :)
    return b.isinside(pv);
  },"points"_a, R"pbdoc(
    Determine whether each of the provided reciprocal lattice points is located
    within the first Brillouin zone

    Parameters
    ----------
    Q : :py:class:`numpy.ndarray`
      A 2 dimensional array of three-vectors (``Q.shape[1]==3``) expressed in
      units of the reciprocal lattice.

    Returns
    -------
    :py:class:`numpy.ndarray`
      One dimensional logical array with ``True`` indicating 'inside'
  )pbdoc");

  cls.def("moveinto",[](CLS &b, py::array_t<double> Q, int threads){
    // brille::Array<double> sp = brille::py2a(Q);
    brille::Array2<double> sp = brille::py2a2(Q);
    if (sp.size(sp.ndim()-1) != 3)
      throw std::runtime_error("The last dimension must have size 3");
    auto Qv = LQVec<double>(b.get_lattice(),  sp); // view
    auto qv = LQVec<double>(b.get_lattice(), sp.shape(), sp.stride()); // output
    auto tauv = LQVec<int>(b.get_lattice(), sp.shape(), sp.stride()); // output
    bool success = b.moveinto(Qv,qv,tauv,threads);
    if (!success) throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
    return py::make_tuple(brille::a2py(qv), brille::a2py(tauv));
  }, "Q"_a, "threads"_a=0, R"pbdoc(
    Find points equivalent to those provided within the first Brillouin zone.

    Parameters
    ----------
    Q : :py:class:`numpy.ndarray`
      A 2 dimensional array of three-vectors (``Q.shape[1]==3``) expressed in
      units of the reciprocal lattice.
    threads : integer, optional
      The number of parallel threads that should be used. If this value is less
      than one the maximum number of OpenMP threads will be used -- this value
      can be controlled by the environment variable ``OMP_NUM_THREADS`` and is
      typically the number of logical cores if not explicitly set.

    Returns
    -------
    :py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`
      The floating point array of equivalent reduced :math:`\mathbf{q}`
      points for all :math:`\mathbf{Q}`, and an integer array filled with
      :math:`\boldsymbol{\tau} = \mathbf{Q}-\mathbf{q}`.
  )pbdoc");

  cls.def("ir_moveinto",[](CLS &b, py::array_t<double> Q, int threads){
    // brille::Array<double> sp = brille::py2a(Q);
    brille::Array2<double> sp = brille::py2a2(Q);
    if (sp.size(sp.ndim()-1) != 3)
      throw std::runtime_error("The last dimension must have size 3");
    auto Qv = LQVec<double>(b.get_lattice(),  sp); // view
    // prepare intermediate outputs
    auto qv = LQVec<double>(b.get_lattice(), sp.shape(), sp.stride()); // output
    auto tauv = LQVec<int>(b.get_lattice(), sp.shape(), sp.stride()); // output
    std::vector<size_t> rotidx(Qv.numel()/3), invrotidx(Qv.numel()/3);
    if (!b.ir_moveinto(Qv, qv, tauv, rotidx, invrotidx, threads))
      throw std::runtime_error("Moving points into irreducible zone failed.");
    // get the pointgroup symmetry operations indexed by rotidx and invrotidx
    PointSymmetry ptsym = b.get_pointgroup_symmetry();
    // prepare Python outputs
    // The rotations array has an extra dimension compared to q and tau
    std::vector<pybind11::ssize_t> sh;
    for (auto s: Qv.shape()) sh.push_back(static_cast<pybind11::ssize_t>(s));
    sh.push_back(3);
    auto rout = py::array_t<int,    py::array::c_style>(sh);
    auto invrout = py::array_t<int, py::array::c_style>(sh);
    // grab pointers to the underlying data blocks
    int *rptr = (int *) rout.request().ptr;
    int *iptr = (int *) invrout.request().ptr;
    for (size_t i=0; i<Qv.numel()/3; ++i)
      for (size_t j=0; j<3u; ++j) for (size_t k=0; k<3u; ++k) {
        rptr[9u*i+3u*j+k] = ptsym.get(rotidx[i])[3u*j+k];
        iptr[9u*i+3u*j+k] = ptsym.get(invrotidx[i])[3u*j+k];
      }
    return py::make_tuple(brille::a2py(qv), brille::a2py(tauv), rout, invrout);
  }, "Q"_a, "threads"_a=0, R"pbdoc(
    Find points equivalent to those provided within the irreducible Brillouin zone.

    The BrillouinZone object defines a volume of reciprocal space which contains
    an irreducible part of the full reciprocal-space. This method will find
    points equivalent under the operations of the lattice which fall within this
    irreducible volume.

    Parameters
    ----------
    Q : :py:class:`numpy.ndarray`
      A 2 dimensional array of three-vectors (``Q.shape[1]==3``) expressed in
      units of the reciprocal lattice.
    threads : integer, optional
      The number of parallel threads that should be used. If this value is less
      than one the maximum number of OpenMP threads will be used -- this value
      can be controlled by the environment variable ``OMP_NUM_THREADS`` and is
      typically the number of logical cores if not explicitly set.

    Returns
    -------
    Qir : :py:class:`numpy.ndarray`
      The array of equivalent irreducible :math:`\mathbf{q}_\text{ir}` points
      for all :math:`\mathbf{Q}`;
    tau : :py:class:`numpy.ndarray`
      the closest reciprocal lattice vector, :math:`\boldsymbol{\tau}`,
      to each :math:`\mathbf{Q}`;
    R : :py:class:`numpy.ndarray`
      the pointgroup symmetry operation :math:`R`
    Rinv : :py:class:`numpy.ndarray`
      the inverse point group symmetry operation which obey
      :math:`\mathbf{Q} = R^{-1} \mathbf{q}_\text{ir} + \boldsymbol{\tau}`.
  )pbdoc");

  cls.def("ir_moveinto_wedge",[](CLS &b, py::array_t<double> Q, int threads){
    // brille::Array<double> sp = brille::py2a(Q);
    brille::Array2<double> sp = brille::py2a2(Q);
    if (sp.size(sp.ndim()-1) != 3)
      throw std::runtime_error("The last dimension must have size 3");
    auto Qv = LQVec<double>(b.get_lattice(), sp); // view
    // prepare intermediate outputs
    auto qv = LQVec<double>(b.get_lattice(), sp.shape(), sp.stride()); // output
    std::vector<std::array<int,9>> rots(Qv.numel()/3);
    std::vector<size_t> ridx(Qv.numel()/3);
    if (!b.ir_moveinto_wedge(Qv, qv, ridx, threads))
      throw std::runtime_error("Moving points into irreducible zone failed.");
    // prepare Python outputs
    // The rotations array has an extra dimension compared to q and tau
    std::vector<pybind11::ssize_t> sh;
    for (auto s: Qv.shape()) sh.push_back(static_cast<pybind11::ssize_t>(s));
    sh.push_back(3);
    auto rout = py::array_t<int,    py::array::c_style>(sh);
    PointSymmetry ptsym = b.get_pointgroup_symmetry();
    // grab pointers to the underlying data blocks
    int *rptr = (int *) rout.request().ptr;
    for (size_t i=0; i<Qv.numel()/3; ++i)
      for (size_t j=0; j<3u; ++j) for (size_t k=0; k<3u; ++k)
        rptr[9u*i+3u*j+k] = ptsym.get(ridx[i])[3u*j+k];
    return py::make_tuple(brille::a2py(qv), rout);
  }, "Q"_a, "threads"_a=0, R"pbdoc(
    Find points equivalent to those provided within the irreducible wedge.

    The BrillouinZone object defines a wedge of reciprocal space which contains
    an irreducible part of the full-space 4π steradian solid angle. This method
    will find points equivalent under the pointgroup operations of the lattice
    which fall within this irreducible solid angle and maintain their absolute
    magnitude.

    Parameters
    ----------
    Q : :py:class:`numpy.ndarray`
      A 2 dimensional array of three-vectors (``Q.shape[1]==3``) expressed in
      units of the reciprocal lattice.
    threads : integer, optional (default 0)
      The number of parallel threads that should be used. If this value is less
      than one the maximum number of OpenMP threads will be used -- this value
      can be controlled by the environment variable ``OMP_NUM_THREADS`` and is
      typically the number of logical cores if not explicitly set.

    Returns
    -------
    :py:class:`numpy.ndarray`, :py:class:`numpy.ndarray`
      The array of equivalent in-wedge :math:`\mathbf{Q}_\text{ir}` points
      for all :math:`\mathbf{Q}`, and the pointgroup operation fulfilling
      :math:`\mathbf{Q}_\text{ir} = R \mathbf{Q}`.
  )pbdoc");

#ifdef USE_HIGHFIVE
  const std::string default_entry("BrillouinZone");
  const std::string default_flags("ac");
  cls.def("to_file",[](CLS& cobj, const std::string& filename, const std::string& entry, const std::string& flags){
        using namespace HighFive;
        unsigned flag{0u};
        if (flags.find('r') != std::string::npos) flag |= File::ReadOnly;
        if (flags.find('x') != std::string::npos) flag |= File::Excl;
        if (flags.find('a') != std::string::npos) flag |= File::ReadWrite;
        if (flags.find('c') != std::string::npos) flag |= File::Create;
        if (flags.find('t') != std::string::npos) flag |= File::Truncate;
        info_update("Provided flags", flags," is translated to ",flag);
        return cobj.to_hdf(filename, entry, flag);
      }, "filename"_a, "entry"_a=default_entry, "flags"_a=default_flags,
      R"pbdoc(
  Save the object to an HDF5 file

  Parameters
  ----------
  filename : str
    The full path specification for the file to write into
  entry: str
    The group path, e.g., "my/cool/bz", where to write inside the file,
    with a default equal to BrillouinZone name
  flags: str
    The HDF5 permissions to use when opening the file. Default 'a' writes to an
    existing file -- if `entry` exists in the file it is overwritten.

  Note
  ----
  Possible `flags` are:

  +---------+-------------------------+----------------+
  | `flags` | meaning                 | HDF equivalent |
  +=========+=========================+================+
  | 'r'     | read                    | H5F_ACC_RDONLY |
  +---------+-------------------------+----------------+
  | 'x'     | write, error if exists  | H5F_ACC_EXCL   |
  +---------+-------------------------+----------------+
  | 'a'     | write, append to file   | H5F_ACC_RDWR   |
  +---------+-------------------------+----------------+
  | 'c'     | write, error if exists  | H5F_ACC_CREAT  |
  +---------+-------------------------+----------------+
  | 't'     | write, replace existing | H5F_ACC_TRUNC  |
  +---------+-------------------------+----------------+


  Returns
  -------
  bool
    Indication of writing success.

  )pbdoc");

  // how do we define this static?
  cls.def_static("from_file",[](const std::string& filename, const std::string& entry){
        return CLS::from_hdf(filename, entry);
      }, "filename"_a, "entry"_a=default_entry,
      R"pbdoc(
  Save the object to an HDF5 file

  Parameters
  ----------
  filename : str
    The full path specification for the file to read from
  entry: str
    The group path, e.g., "my/cool/bz", where to read from inside the file,
    with a default equal to the object Class name

  Returns
  -------
  clsObj

  )pbdoc");
#endif //USE_HIGHFIVE
}
