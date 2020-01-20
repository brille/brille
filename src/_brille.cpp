/* Copyright 2019 Greg Tucker
//
// This file is part of brille.
//
// brille is free software: you can redistribute it and/or modify it under the
// terms of the GNU Affero General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// brille is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with brille. If not, see <https://www.gnu.org/licenses/>.            */


#include "version_info.hpp"
#include "pointgroup.hpp"
#include "_c_to_python.hpp"
#include "_grid.hpp"
#include "_lattice.hpp"
#include "_mesh.hpp"
#include "_polyhedron.hpp"
#include "_trellis.hpp"
#include "_nest.hpp"
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "hall_symbol.hpp"

namespace py = pybind11;
using namespace pybind11::literals; // bring in "[name]"_a to be interpreted as py::arg("[name]")

typedef long slong; // ssize_t is only defined for gcc?

PYBIND11_MODULE(_brille,m){
  m.doc() = "brille for dealing with symmetry of the first irreducible Brillouin Zone";

  {
  using namespace brille::version;
  m.attr("__version__") = version_number;
  m.attr("version") = long_version();
  m.attr("git_revision") = git_revision;
  m.attr("build_datetime") = build_datetime;
  m.attr("build_hostname") = build_hostname;
  }

  py::class_<Lattice>(m,"Lattice")
    .def(py::init<double,double,double,double,double,double,int>(),
         py::arg("a"),py::arg("b"),py::arg("c"),
         py::arg("alpha")=PI/2,py::arg("beta")=PI/2,py::arg("gamma")=PI/2,
         py::arg("HallNumber")=1)
    .def(py::init([](py::array_t<double> lens, py::array_t<double> angs, int hall) {
      py::buffer_info linfo = lens.request(), ainfo = angs.request();
      if ( linfo.ndim!=1 || ainfo.ndim!=1)
        throw std::runtime_error("Number of dimensions must be one");
      if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
        throw std::runtime_error("(At least) three lengths and angles required.");
      return Lattice((double*) linfo.ptr, linfo.strides,
                     (double*) ainfo.ptr, ainfo.strides, hall);
    }),py::arg("lengths"),py::arg("angles"),py::arg("HallNumber")=1)
    .def(py::init([](py::array_t<double> vecs, int hall) {
      py::buffer_info info = vecs.request();
      if ( info.ndim!=2 )
        throw std::runtime_error("Number of dimensions must be two");
      if ( info.shape[0] != 3 || info.shape[1] != 3 )
        throw std::runtime_error("Three three-vectors required.");
      return Lattice((double*) info.ptr, info.strides, hall);
    }),py::arg("vectors"),py::arg("HallNumber")=1)
    .def_property_readonly("a",     &Lattice::get_a)
    .def_property_readonly("b",     &Lattice::get_b)
    .def_property_readonly("c",     &Lattice::get_c)
    .def_property_readonly("alpha", &Lattice::get_alpha)
    .def_property_readonly("beta",  &Lattice::get_beta)
    .def_property_readonly("gamma", &Lattice::get_gamma)
    .def_property_readonly("volume",&Lattice::get_volume)
    .def_property("hall",&Lattice::get_hall,&Lattice::set_hall)
    .def_property_readonly("spacegroup",&Lattice::get_spacegroup_object)
    .def("get_covariant_metric_tensor",[](Lattice &l){
      auto result = py::array_t<double, py::array::c_style >({3,3});
      py::buffer_info bi = result.request();
      double *cmt = (double *) bi.ptr;
      l.get_covariant_metric_tensor( cmt );
      return result;
    })
    .def("get_contravariant_metric_tensor",[](Lattice &l){
      auto result = py::array_t<double, py::array::c_style >({3,3});
      py::buffer_info bi = result.request();
      l.get_contravariant_metric_tensor((double *) bi.ptr );
      return result;
    })
    .def_property_readonly("star",[](Lattice &){throw std::runtime_error("Bare Lattices do not have a reciprocal!");})
    .def("issame",&Lattice::issame)
    .def("__eq__",&Lattice::issame)
    .def("isapprox",&Lattice::isapprox)
    .def("isstar",[](Lattice &, Lattice ){throw std::runtime_error("Bare Lattices do not have a reciprocal!");})
    .def("__repr__",&Lattice::string_repr)
    ;

  py::class_<Direct,Lattice> direct(m,"Direct");
  py::class_<Reciprocal,Lattice> reciprocal(m,"Reciprocal");

  declare_lattice_methods<Direct>(direct,"Å");
  declare_lattice_methods<Reciprocal>(reciprocal,"Å⁻¹");
  reciprocal.def("get_B_matrix",[](Reciprocal &r){
    auto result = py::array_t<double, py::array::c_style>({3,3});
    py::buffer_info bi = result.request();
    r.get_B_matrix((double *)bi.ptr);
    return result;
  });
  reciprocal.def_property_readonly("B",[](Reciprocal &r){
    auto result = py::array_t<double, py::array::c_style>({3,3});
    py::buffer_info bi = result.request();
    r.get_B_matrix((double *)bi.ptr);
    return result;
  });

  py::class_<BrillouinZone> bz(m,"BrillouinZone");
  bz.def(py::init<Reciprocal,bool,int,bool,bool>(),
         py::arg("lattice"),
         py::arg("use_primitive")=true,
         py::arg("search_length")=1,
         py::arg("time_reversal_symmetry")=false,
         py::arg("wedge_search")=true
       )
    .def_property_readonly("lattice", [](const BrillouinZone &b){ return b.get_lattice();} )
    // access the polyhedra directly
    .def_property_readonly("polyhedron",&BrillouinZone::get_polyhedron)
    .def_property_readonly("ir_polyhedron",[](const BrillouinZone &b){return b.get_ir_polyhedron(true);})
    .def_property_readonly("ir_polyhedron_generated",[](const BrillouinZone &b){return b.get_ir_polyhedron(false);})
    // first Brillouin zone polyhedron
    .def_property_readonly("normals",                  [](const BrillouinZone &b){return av2np(b.get_normals().get_hkl());})
    .def_property_readonly("normals_invA",             [](const BrillouinZone &b){return av2np(b.get_normals().get_xyz());})
    .def_property_readonly("normals_primitive",        [](const BrillouinZone &b){return av2np(b.get_primitive_normals().get_hkl());})
    .def_property_readonly("points",                   [](const BrillouinZone &b){return av2np(b.get_points().get_hkl());})
    .def_property_readonly("points_invA",              [](const BrillouinZone &b){return av2np(b.get_points().get_xyz());})
    .def_property_readonly("points_primitive",         [](const BrillouinZone &b){return av2np(b.get_primitive_points().get_hkl());})
    .def_property_readonly("vertices",                 [](const BrillouinZone &b){return av2np(b.get_vertices().get_hkl());})
    .def_property_readonly("vertices_invA",            [](const BrillouinZone &b){return av2np(b.get_vertices().get_xyz());})
    .def_property_readonly("half_edge_points",         [](const BrillouinZone &b){return av2np(b.get_half_edges().get_hkl());})
    .def_property_readonly("half_edge_points_invA",    [](const BrillouinZone &b){return av2np(b.get_half_edges().get_xyz());})
    .def_property_readonly("vertices_primitive",       [](const BrillouinZone &b){return av2np(b.get_primitive_vertices().get_hkl());})
    .def_property_readonly("faces_per_vertex",  &BrillouinZone::get_faces_per_vertex)
    .def_property_readonly("vertices_per_face", &BrillouinZone::get_vertices_per_face)
    // irreducible first Brillouin zone polyhedron
    .def_property_readonly("ir_normals",                  [](const BrillouinZone &b){return av2np(b.get_ir_normals().get_hkl());})
    .def_property_readonly("ir_normals_invA",             [](const BrillouinZone &b){return av2np(b.get_ir_normals().get_xyz());})
    .def_property_readonly("ir_normals_primitive",        [](const BrillouinZone &b){return av2np(b.get_ir_primitive_normals().get_hkl());})
    .def_property_readonly("ir_points",                   [](const BrillouinZone &b){return av2np(b.get_ir_points().get_hkl());})
    .def_property_readonly("ir_points_invA",              [](const BrillouinZone &b){return av2np(b.get_ir_points().get_xyz());})
    .def_property_readonly("ir_points_primitive",         [](const BrillouinZone &b){return av2np(b.get_ir_primitive_points().get_hkl());})
    .def_property_readonly("ir_vertices",                 [](const BrillouinZone &b){return av2np(b.get_ir_vertices().get_hkl());})
    .def_property_readonly("ir_vertices_invA",            [](const BrillouinZone &b){return av2np(b.get_ir_vertices().get_xyz());})
    .def_property_readonly("ir_vertices_primitive",       [](const BrillouinZone &b){return av2np(b.get_ir_primitive_vertices().get_hkl());})
    .def_property_readonly("ir_faces_per_vertex",  &BrillouinZone::get_ir_faces_per_vertex)
    .def_property_readonly("ir_vertices_per_face", &BrillouinZone::get_ir_vertices_per_face)
    // irreducible reciprocal space wedge
    .def_property_readonly("wedge_normals",            [](const BrillouinZone &b){return av2np(b.get_ir_wedge_normals().get_hkl());})
    .def_property_readonly("wedge_normals_invA",       [](const BrillouinZone &b){return av2np(b.get_ir_wedge_normals().get_xyz());})
    .def_property_readonly("wedge_normals_primitive",  [](const BrillouinZone &b){return av2np(b.get_primitive_ir_wedge_normals().get_hkl());})
    .def("isinside",[](BrillouinZone &b, py::array_t<double> p){
      py::buffer_info bi = p.request();
      ssize_t ndim = bi.ndim;
      if ( bi.shape[ndim-1] !=3) throw std::runtime_error("one or more 3-dimensional points is required");
      ssize_t npts = 1;
      if (ndim > 1)  for (ssize_t i=0; i<ndim-1; i++) npts *= bi.shape[i];
      LQVec<double> pv( b.get_lattice(), (double*) bi.ptr, bi.shape, bi.strides); // this is a copy :(
      ArrayVector<bool> resultv = b.isinside(pv);
      std::vector<ssize_t> outshape(ndim>1 ? ndim-1 : 1);
      if (ndim > 1){
        for (ssize_t i=0; i<ndim-1; i++) outshape[i] = bi.shape[i];
      } else {
        outshape[0] = 1;
      }
      auto result = py::array_t<bool, py::array::c_style>(outshape);
      bool *rptr = (bool *) result.request().ptr;
      for (ssize_t i=0; i<npts; i++) rptr[i] = resultv.getvalue(i);
      return result;
    },py::arg("points"))
    .def("moveinto",[](BrillouinZone &b, py::array_t<double> Q){
      py::buffer_info bi = Q.request();
      ssize_t ndim=bi.ndim;
      if (bi.shape[ndim-1] !=3) throw std::runtime_error("one or more 3-dimensional Q points is required");
      ssize_t npts = 1;
      if (ndim > 1) for (ssize_t i=0; i<ndim-1; i++) npts *= bi.shape[i];
      LQVec<double> Qv( b.get_lattice(),(double*) bi.ptr, bi.shape, bi.strides); // memcopy
      LQVec<double> qv(b.get_lattice(), npts); // output
      LQVec<int>  tauv(b.get_lattice(), npts); // output
      bool success = b.moveinto(Qv,qv,tauv);
      if (!success) throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
      auto qout = py::array_t<double, py::array::c_style>(bi.shape);
      auto tout = py::array_t<int, py::array::c_style>(bi.shape);
      double *qptr = (double *) qout.request().ptr;
      int *tptr = (int *) tout.request().ptr;
      for (ssize_t i=0; i<npts; ++i){
        for (size_t j=0; j<3u; ++j){
          qptr[3u*i+j] = qv.getvalue(i,j);
          tptr[3u*i+j] = tauv.getvalue(i,j);
        }
      }
      return py::make_tuple(qout,tout);
    }, py::arg("Q"))
    .def("ir_moveinto",[](BrillouinZone &b, py::array_t<double> Q){
      py::buffer_info bi = Q.request();
      ssize_t ndim = bi.ndim;
      if (bi.shape[ndim-1] != 3)
        throw std::runtime_error("One or more 3-dimensional Q points are required.");
      ssize_t npts = 1;
      if (ndim > 1) for (ssize_t i=0; i<ndim-1; ++i) npts *= bi.shape[i];
      LQVec<double> Qv(b.get_lattice(), (double*) bi.ptr, bi.shape, bi.strides);
      // prepare intermediate outputs
      LQVec<double> qv(b.get_lattice(), npts);
      LQVec<int>  tauv(b.get_lattice(), npts);
      std::vector<std::array<int,9>> rots(npts);
      if (!b.ir_moveinto(Qv, qv, tauv, rots))
        throw std::runtime_error("Moving points into irreducible zone failed.");
      // prepare Python outputs
      auto qout = py::array_t<double, py::array::c_style>(bi.shape);
      auto tout = py::array_t<int,    py::array::c_style>(bi.shape);
      // The rotations array has an extra dimension compared to q and tau
      bi.shape.push_back(3);
      auto rout = py::array_t<int,    py::array::c_style>(bi.shape);
      // grab pointers to the underlying data blocks
      double *qptr = (double *) qout.request().ptr;
      int *tptr = (int *) tout.request().ptr;
      int *rptr = (int *) rout.request().ptr;
      for (ssize_t i=0; i<npts; ++i)
      for (size_t j=0; j<3u; ++j){
        qptr[3u*i+j] = qv.getvalue(i,j);
        tptr[3u*i+j] = tauv.getvalue(i,j);
        for (size_t k=0; k<3u; ++k) rptr[9u*i+3u*j+k] = rots[i][3u*j+k];
      }
      return py::make_tuple(qout, tout, rout);
    }, py::arg("Q"))
    .def("ir_moveinto_wedge",[](BrillouinZone &b, py::array_t<double> Q){
      py::buffer_info bi = Q.request();
      ssize_t ndim = bi.ndim;
      if (bi.shape[ndim-1] != 3)
        throw std::runtime_error("One or more 3-dimensional Q points are required.");
      ssize_t npts = 1;
      if (ndim > 1) for (ssize_t i=0; i<ndim-1; ++i) npts *= bi.shape[i];
      LQVec<double> Qv(b.get_lattice(), (double*) bi.ptr, bi.shape, bi.strides);
      // prepare intermediate outputs
      LQVec<double> qv(b.get_lattice(), npts);
      std::vector<std::array<int,9>> rots(npts);
      if (!b.ir_moveinto_wedge(Qv, qv, rots))
        throw std::runtime_error("Moving points into irreducible zone failed.");
      // prepare Python outputs
      auto qout = py::array_t<double, py::array::c_style>(bi.shape);
      // The rotations array has an extra dimension compared to q and tau
      bi.shape.push_back(3);
      auto rout = py::array_t<int,    py::array::c_style>(bi.shape);
      // grab pointers to the underlying data blocks
      double *qptr = (double *) qout.request().ptr;
      int *rptr = (int *) rout.request().ptr;
      for (ssize_t i=0; i<npts; ++i)
      for (size_t j=0; j<3u; ++j){
        qptr[3u*i+j] = qv.getvalue(i,j);
        for (size_t k=0; k<3u; ++k) rptr[9u*i+3u*j+k] = rots[i][3u*j+k];
      }
      return py::make_tuple(qout, rout);
    }, py::arg("Q"));

    declare_bzgridq<double>(m,"");
    declare_bzgridqe<double>(m,"");
    declare_bzgridq<std::complex<double>>(m,"complex");
    declare_bzgridqe<std::complex<double>>(m,"complex");

    declare_bzmeshq<double>(m,"");
    declare_bzmeshq<std::complex<double>>(m,"complex");

    declare_bztrellisq<double>(m,"");
    declare_bztrellisq<std::complex<double>>(m,"complex");

    declare_bznestq<double>(m,"");
    declare_bznestq<std::complex<double>>(m,"complex");

    py::class_<PrimitiveTransform> pt(m,"PrimitiveTransform");
    pt.def(py::init<int>(),py::arg("Hall number"));
    pt.def_property_readonly("P",[](const PrimitiveTransform &p){
      std::vector<ssize_t> sz={3,3};
      return sa2np(sz,p.get_P());
    });
    pt.def_property_readonly("invP",[](const PrimitiveTransform &p){
      std::vector<ssize_t> sz={3,3};
      return sa2np(sz,p.get_invP());
    });
    pt.def_property_readonly("Pt",[](const PrimitiveTransform &p){
      std::vector<ssize_t> sz={3,3};
      return sa2np(sz,p.get_Pt());
    });
    pt.def_property_readonly("invPt",[](const PrimitiveTransform &p){
      std::vector<ssize_t> sz={3,3};
      return sa2np(sz,p.get_invPt());
    });
    pt.def_property_readonly("does_anything",&PrimitiveTransform::does_anything);
    pt.def_property_readonly("is_primitive",&PrimitiveTransform::does_nothing);
    pt.def("__repr__",&PrimitiveTransform::string_repr);

    py::class_<Spacegroup> spg(m,"Spacegroup");
    spg.def(py::init<int>(),py::arg("Hall number"));
    spg.def_property_readonly("hall_number", &Spacegroup::get_hall_number);
    spg.def_property_readonly("international_table_number", &Spacegroup::get_international_table_number);
    spg.def_property_readonly("pointgroup_number", &Spacegroup::get_pointgroup_number);
    spg.def_property_readonly("schoenflies_symbol", &Spacegroup::get_schoenflies_symbol);
    spg.def_property_readonly("hall_symbol", &Spacegroup::get_hall_symbol);
    spg.def_property_readonly("international_table_symbol", &Spacegroup::get_international_table_symbol);
    spg.def_property_readonly("international_table_full", &Spacegroup::get_international_table_full);
    spg.def_property_readonly("international_table_short", &Spacegroup::get_international_table_short);
    spg.def_property_readonly("choice", &Spacegroup::get_choice);
    spg.def("__repr__",&Spacegroup::string_repr);

    py::class_<Pointgroup> pg(m, "Pointgroup");
    pg.def(py::init<int>(), py::arg("Pointgroup number"));
    pg.def_property_readonly("number",&Pointgroup::get_number);
    pg.def_property_readonly("symbol",&Pointgroup::get_symbol);
    pg.def_property_readonly("holohedry",&Pointgroup::get_holohedry_string);
    pg.def_property_readonly("laue",&Pointgroup::get_laue_string);
    pg.def("__repr__",&Pointgroup::to_string);

    py::class_<Symmetry> sym(m, "Symmetry");
    sym.def(py::init([](int hall){return Spacegroup(hall).get_spacegroup_symmetry();}),py::arg("Hall number"));
    sym.def_property_readonly("size",&Symmetry::size);
    sym.def_property_readonly("W",[](Symmetry& ps){
      std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3, 3};
      return sva2np(sz, ps.getallr());
    });
    sym.def_property_readonly("w",[](Symmetry& ps){
      std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3};
      return sva2np(sz, ps.getallt());
    });
    sym.def("generate",&Symmetry::generate);
    sym.def("generators",&Symmetry::generators);
    sym.def(py::self == py::self);

    py::class_<PointSymmetry> psym(m, "PointSymmetry");
    psym.def(py::init( [](int hall, int time_reversal){return Spacegroup(hall).get_pointgroup_symmetry(time_reversal);}),py::arg("Hall_number"),py::arg("time_reversal")=0);
    psym.def_property_readonly("size",&PointSymmetry::size);
    psym.def_property_readonly("W",[](PointSymmetry& ps){
      std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3, 3};
      return sva2np(sz, ps.getall());
    });
    psym.def_property_readonly("order",[](PointSymmetry& ps){
      std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size())};
      return sv2np(sz, ps.orders());
    });
    psym.def_property_readonly("isometry",[](PointSymmetry& ps){
      std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size())};
      return sv2np(sz, ps.isometries());
    });
    psym.def_property_readonly("axis",[](PointSymmetry& ps){
      std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3};
      return sva2np(sz, ps.axes());
    });
    psym.def_property_readonly("generate",&PointSymmetry::generate);
    psym.def_property_readonly("generators",&PointSymmetry::generators);
    psym.def("nfolds",&PointSymmetry::nfolds);

    declare_polyhedron<Polyhedron>(m, "");


    py::class_<HallSymbol> hsbl(m, "HallSymbol");
    hsbl.def(py::init([](std::string symbol){
      HallSymbol out;
      out.from_ascii(symbol);
      return out;
    }),py::arg("Hall_symbol"));
    hsbl.def("__repr__",&HallSymbol::to_ascii);
    hsbl.def_property_readonly("generators",&HallSymbol::get_generators);
}
