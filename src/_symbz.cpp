#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <vector>
#include <complex>
#include <omp.h>

// #include <symbz.h>
#include "symbz.h"
#include "arithmetic.h"
#include "linear_algebra.h"
#include "lattice.h"
#include "pointgroup.h"
#include "spg_database.h"
#include "symmetry.h"
#include "bz.h"
#include "grid.h"
#include "bz_grid.h"

#include "version_info.h"

#include "_binding.h"

namespace py = pybind11;
using namespace pybind11::literals; // bring in "[name]"_a to be interpreted as py::arg("[name]")

typedef long slong; // ssize_t is only defined for gcc?

PYBIND11_MODULE(_symbz,m){
	m.doc() = "SymBZ for dealing with symmetry of the first Brillouin Zone";

	{
	using namespace symbz::version;
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
		.def(py::init([](py::array_t<double> lens, py::array_t<double> angs) {
			py::buffer_info linfo = lens.request(), ainfo = angs.request();
			if ( linfo.ndim!=1 || ainfo.ndim!=1)
				throw std::runtime_error("Number of dimensions must be one");
			if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
				throw std::runtime_error("(At least) three lengths and angles required.");
			double *lengths = (double *) linfo.ptr, *angles = (double *) ainfo.ptr;
			return Lattice(lengths,angles);
		}))
		.def_property_readonly("a",     &Lattice::get_a)
		.def_property_readonly("b",     &Lattice::get_b)
		.def_property_readonly("c",     &Lattice::get_c)
		.def_property_readonly("alpha", &Lattice::get_alpha)
		.def_property_readonly("beta",  &Lattice::get_beta)
		.def_property_readonly("gamma", &Lattice::get_gamma)
		.def_property_readonly("volume",&Lattice::get_volume)
		.def_property("hall",&Lattice::get_hall,&Lattice::set_hall)
		.def("fill_covariant_metric_tensor",[](Lattice &l, py::array_t<double> cmt){
			py::buffer_info bi = cmt.request();
			if (bi.ndim!=2) throw std::runtime_error("Number of dimensions must be 2");
			if (bi.shape[0] !=3 || bi.shape[1] != 3) throw std::runtime_error("Array must be 3x3");
			l.get_covariant_metric_tensor( (double *) bi.ptr);
		})
		.def("get_covariant_metric_tensor",[](Lattice &l){
			auto result = py::array_t<double, py::array::c_style >({3,3});
			py::buffer_info bi = result.request();
			double *cmt = (double *) bi.ptr;
			l.get_covariant_metric_tensor( cmt );
			return result;
		})
		.def("fill_contravariant_metric_tensor",[](Lattice &l, py::array_t<double> cmt){
			py::buffer_info bi = cmt.request();
			if (bi.ndim!=2) throw std::runtime_error("Number of dimensions must be 2");
			if (bi.shape[0] !=3 || bi.shape[1] != 3) throw std::runtime_error("Array must be 3x3");
			l.get_contravariant_metric_tensor( (double *) bi.ptr);
		})
		.def("get_contravariant_metric_tensor",[](Lattice &l){
			auto result = py::array_t<double, py::array::c_style >({3,3});
			py::buffer_info bi = result.request();
			l.get_contravariant_metric_tensor((double *) bi.ptr );
			return result;
		})
		.def("star",[](Lattice &l){throw std::runtime_error("Bare Lattices do not have a reciprocal!");})
		.def("issame",&Lattice::issame)
		.def("__eq__",&Lattice::issame)
		.def("isstar",[](Lattice &l, Lattice a){throw std::runtime_error("Bare Lattices do not have a reciprocal!");})
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

	py::class_<BrillouinZone> bz(m,"BrillouinZone");
	bz.def(py::init<Reciprocal,int>(),py::arg("lattice"),py::arg("search_length")=1)
	  .def_property_readonly("lattice", [](const BrillouinZone &b){ return b.get_lattice();} )
	  .def_property_readonly("vertices",    [](const BrillouinZone &b){return av2np(b.get_vertices().get_hkl());})
		.def_property_readonly("vertices_invA",[](const BrillouinZone &b){return av2np(b.get_vertices().get_xyz());})
		.def_property_readonly("faces",    [](const BrillouinZone &b){return av2np(b.get_faces().get_hkl());})
		.def_property_readonly("faces_invA",[](const BrillouinZone &b){return av2np(b.get_faces().get_xyz());})
		.def_property_readonly("faces_per_vertex",    [](const BrillouinZone &b){return av2np(b.get_faces_per_vertex());})
		.def("isinside",[](BrillouinZone &b, py::array_t<double, py::array::c_style> p){
			py::buffer_info bi = p.request();
			ssize_t ndim = bi.ndim;
			if ( bi.shape[ndim-1] !=3) throw std::runtime_error("one or more 3-dimensional points is required");
			ssize_t npts = 1;
			if (ndim > 1)	for (ssize_t i=0; i<ndim-1; i++) npts *= bi.shape[i];
			LQVec<double> pv( b.get_lattice(), npts, (double*) bi.ptr ); // this is a copy :(
			ArrayVector<bool> resultv = b.isinside(pv);
			std::vector<ssize_t> outshape(ndim>1 ? ndim-1 : 1);
			if (ndim > 1){
				for (ssize_t i=0; i<ndim-1; i++) outshape[i] = bi.shape[i];
			} else {
				outshape[0] = 1;
			}
			auto result = py::array_t<bool, py::array::c_style>(outshape);
			bool *rptr = (bool *) result.request().ptr;
			for (size_t i=0; i<npts; i++) rptr[i] = resultv.getvalue(i);
			return result;
		},py::arg("points"))
		.def("moveinto",[](BrillouinZone &b, py::array_t<double, py::array::c_style> Q){
			py::buffer_info bi = Q.request();
			ssize_t ndim=bi.ndim;
			if (bi.shape[ndim-1] !=3) throw std::runtime_error("one or more 3-dimensional Q points is required");
			ssize_t npts = 1;
			if (ndim > 1) for (ssize_t i=0; i<ndim-1; i++) npts *= bi.shape[i];
			LQVec<double> Qv( b.get_lattice(), npts, (double*) bi.ptr); // memcopy
			LQVec<double> qv(b.get_lattice(), npts); // output
			LQVec<int>  tauv(b.get_lattice(), npts); // output
			bool success = b.moveinto(Qv,qv,tauv);
			if (!success) throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
			auto qout = py::array_t<double, py::array::c_style>(bi.shape);
			auto tout = py::array_t<int, py::array::c_style>(bi.shape);
			double *qptr = (double *) qout.request().ptr;
			int *tptr = (int *) tout.request().ptr;
			for (size_t i=0; i<npts; ++i){
				for (size_t j=0; j<3u; ++j){
					qptr[3u*i+j] = qv.getvalue(i,j);
					tptr[3u*i+j] = tauv.getvalue(i,j);
				}
			}
			return py::make_tuple(qout,tout);
		}, py::arg("Q"));

		declare_bzgridq<double>(m,"");
		declare_bzgridqe<double>(m,"");
		declare_bzgridq<std::complex<double>>(m,"complex");
		declare_bzgridqe<std::complex<double>>(m,"complex");

}
