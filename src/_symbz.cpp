#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>

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

namespace py = pybind11;
using namespace pybind11::literals; // bring in "[name]"_a to be interpreted as py::arg("[name]")

template<typename T> py::array_t<T> av2np(const ArrayVector<T> av){
	std::vector<ssize_t> shape(2); // ArrayVectors are 2D by default
	shape[0] = av.size();
	shape[1] = av.numel();
	auto np = py::array_t<T,py::array::c_style>(shape);
	T *rptr = (T*) np.request().ptr;
	for (size_t i =0; i< av.size(); i++)
		for (size_t j=0; j< av.numel(); j++)
			rptr[i*av.numel()+j] = av.getvalue(i,j);
	return np;
}


PYBIND11_MODULE(symbz,m){
	m.doc() = "SymBZ for dealing with symmetry of the first Brillouin Zone";

	py::class_<Lattice>(m,"Lattice")
		.def(py::init<double,double,double,double,double,double,double>(),
		     py::arg("a"),py::arg("b"),py::arg("c"),py::arg("alpha"),py::arg("beta"),py::arg("gamma"),py::arg("volume"))
		.def(py::init<double,double,double,double,double,double>(),
		     py::arg("a"),py::arg("b"),py::arg("c"),py::arg("alpha")=PI/2,py::arg("beta")=PI/2,py::arg("gamma")=PI/2)
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

	direct.def("star", &Direct::star, "A real space lattice")
				.def(py::init<double,double,double,double,double,double,double>(),
						 py::arg("a"),py::arg("b"),py::arg("c"),py::arg("alpha"),py::arg("beta"),py::arg("gamma"),py::arg("volume"))
				.def(py::init<double,double,double,double,double,double>(),
						 py::arg("a"),py::arg("b"),py::arg("c"),py::arg("alpha")=PI/2,py::arg("beta")=PI/2,py::arg("gamma")=PI/2)
				.def(py::init( [](py::array_t<double> lens, py::array_t<double> angs) {
					py::buffer_info linfo = lens.request(), ainfo = angs.request();
					if ( linfo.ndim!=1 || ainfo.ndim!=1)
						throw std::runtime_error("Number of dimensions must be one");
					if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
						throw std::runtime_error("(At least) three lengths and angles required.");
					double *lengths = (double *) linfo.ptr, *angles = (double *) ainfo.ptr;
					return Direct(lengths,angles);
				}), py::arg("lengths/Å"), py::arg("angles/radian"))
		    .def("get_xyz_transform",[](Direct &d){
					auto result = py::array_t<double, py::array::c_style>({3,3});
					py::buffer_info bi = result.request();
					d.get_xyz_transform((double *)bi.ptr);
					return result;
				})
				.def("isstar",(bool (Direct::*)(const Direct) const) &Direct::isstar)
				.def("isstar",(bool (Direct::*)(const Reciprocal) const) &Direct::isstar);

	reciprocal.def("star", &Reciprocal::star, "A reciprocal space lattice")
						.def(py::init<double,double,double,double,double,double,double>(),
								 py::arg("a"),py::arg("b"),py::arg("c"),py::arg("alpha"),py::arg("beta"),py::arg("gamma"),py::arg("volume"))
						.def(py::init<double,double,double,double,double,double>(),
								 py::arg("a"),py::arg("b"),py::arg("c"),py::arg("alpha")=PI/2,py::arg("beta")=PI/2,py::arg("gamma")=PI/2)
						.def(py::init([](py::array_t<double> lens, py::array_t<double> angs) {
							py::buffer_info linfo = lens.request(), ainfo = angs.request();
							if ( linfo.ndim!=1 || ainfo.ndim!=1)
								throw std::runtime_error("Number of dimensions must be one");
							if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
								throw std::runtime_error("(At least) three lengths and angles required.");
							double *lengths = (double *) linfo.ptr, *angles = (double *) ainfo.ptr;
							return Reciprocal(lengths,angles);
						}), py::arg("lengths/Å⁻¹"), py::arg("angles/radian"))
						.def("get_B_matrix",[](Reciprocal &r){
							auto result = py::array_t<double, py::array::c_style>({3,3});
							py::buffer_info bi = result.request();
							r.get_B_matrix((double *)bi.ptr);
							return result;
						})
						.def("get_xyz_transform",[](Reciprocal &r){
							auto result = py::array_t<double, py::array::c_style>({3,3});
							py::buffer_info bi = result.request();
							r.get_xyz_transform((double *)bi.ptr);
							return result;
						})
						.def("isstar",(bool (Reciprocal::*)(const Reciprocal) const) &Reciprocal::isstar)
						.def("isstar",(bool (Reciprocal::*)(const Direct) const) &Reciprocal::isstar);

	py::class_<BrillouinZone> bz(m,"BrillouinZone");
	bz.def(py::init<Reciprocal,int>(),py::arg("lattice"),py::arg("search_length")=1)
	  .def_property_readonly("vertices",    [](BrillouinZone &b){return av2np(b.get_vertices().get_hkl());})
		.def_property_readonly("vertices_xyz",[](BrillouinZone &b){return av2np(b.get_vertices().get_xyz());})
		.def_property_readonly("faces",    [](BrillouinZone &b){return av2np(b.get_faces().get_hkl());})
		.def_property_readonly("faces_xyz",[](BrillouinZone &b){return av2np(b.get_faces().get_xyz());})
		.def_property_readonly("faces_per_vertex",    [](BrillouinZone &b){return av2np(b.get_faces());})
		.def("isinside",[](BrillouinZone &b, py::array_t<double, py::array::c_style> p, double tol){
			py::buffer_info bi = p.request();
			ssize_t ndim = bi.ndim;
			if ( bi.shape[ndim-1] !=3) throw std::runtime_error("one or more 3-dimensional points is required");
			ssize_t npts = 1;
			if (ndim > 1)	for (ssize_t i=0; i<ndim-1; i++) npts *= bi.shape[i];
			LQVec<double> pv( b.get_lattice(), npts, (double*) bi.ptr ); // this is a copy :(
			ArrayVector<bool> resultv = b.isinside(&pv,tol);
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
		},py::arg("points"),py::arg("tol")=1e-14)
		.def("moveinto",[](BrillouinZone &b, py::array_t<double, py::array::c_style> Q){
			py::buffer_info bi = Q.request();
			ssize_t ndim=bi.ndim;
			if (bi.shape[ndim-1] !=3) throw std::runtime_error("one or more 3-dimensional Q points is required");
			ssize_t npts = 1;
			if (ndim > 1) for (ssize_t i=0; i<ndim-1; i++) npts *= bi.shape[i];
			LQVec<double> Qv( b.get_lattice(), npts, (double*) bi.ptr); // memcopy
			LQVec<double> qv(b.get_lattice(), npts); // output
			LQVec<int>  tauv(b.get_lattice(), npts); // output
			bool success = b.moveinto(&Qv,&qv,&tauv);
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
		}, py::arg("Q"))
		;

	py::class_<BrillouinZoneGrid3> bzg(m,"BZGrid");
	// Initializer (BrillouinZone, [half-]Number_of_steps vector)
	bzg.def(py::init([](BrillouinZone &b, py::array_t<size_t,py::array::c_style> pyN){
		py::buffer_info bi = pyN.request();
		if (bi.ndim != 1) throw std::runtime_error("N must be a 1-D array");
		if (bi.shape[0] < 3) throw std::runtime_error("N must have three elements");
		BrillouinZoneGrid3 bzg3( b, (size_t*)bi.ptr );
		return bzg3;
	}),py::arg("brillouinzone"),py::arg("N"));
	// Initializer (BrillouinZone, step_size vector, flag_for_whether_step_size_is_in_rlu_or_inverse_angstrom)
  bzg.def(py::init([](BrillouinZone &b, py::array_t<double,py::array::c_style> pyD, const bool& isrlu){
		py::buffer_info bi = pyD.request();
		if (bi.ndim != 1) throw std::runtime_error("stepsize must be a 1-D array");
		if (bi.shape[0] < 3) throw std::runtime_error("stepsize must have three elements");
		return BrillouinZoneGrid3( b, (double*)bi.ptr, isrlu ? 1 : 0 );
	}),py::arg("brillouinzone"),py::arg("step"),py::arg("rlu")=true);
	bzg.def_property_readonly("brillouinzone",[](const BrillouinZoneGrid3& bzg3){ return bzg3.get_brillouinzone();} )
	   .def_property_readonly("hkl",[](const BrillouinZoneGrid3& bzg3){ return av2np(bzg3.get_grid_hkl());} )
		 .def_property_readonly("xyz",[](const BrillouinZoneGrid3& bzg3){ return av2np(bzg3.get_grid_xyz());} )
		 .def_property_readonly("mapped_hkl",[](const BrillouinZoneGrid3& bzg3){ return av2np(bzg3.get_mapped_hkl());} )
		 .def_property_readonly("mapped_xyz",[](const BrillouinZoneGrid3& bzg3){ return av2np(bzg3.get_mapped_xyz());} );
  bzg.def("fill",[](BrillouinZoneGrid3& bzg3, py::array_t<double,py::array::c_style> pydata){
		py::buffer_info bi = pydata.request();
		ssize_t ndim = bi.ndim;
		/* ndim  assumed interpolation data type  numel /
		/   1              scalar                   1   /
		/   2              vector                shape[1]  /
		/   3       matrix / rank 2 tensor       shape[1]*shape[2]       /
		/   N         rank N-1 tensor            prod(shape,1,N-1)      */
		size_t numel=1, numarrays=bi.shape[0];
		if (ndim > 1) for (ssize_t i=1; i<ndim; ++i) numel *= bi.shape[i];
		ArrayVector<double> data(numel, numarrays, (double*)bi.ptr);
		ArrayVector<size_t> shape(1,ndim);
		for (ssize_t i=0; i<ndim; ++i) shape.insert(bi.shape[i], (size_t)i );
		int mapExceedsNewData = bzg3.check_map(data);
		if (mapExceedsNewData) throw std::runtime_error("There are less provided data arrays than unique integers in the mapping.");
		bzg3.replace_data(data,shape); // no error, so this will work for sure
		// return mapExceedsNewData; // let the calling function deal with this?
	});
	bzg.def_property("map",
		/*get map*/ [](BrillouinZoneGrid3& bzg3){
			std::vector<ssize_t> shape(3); // the map is 3D
			shape[0] = bzg3.size(0);
			shape[1] = bzg3.size(1);
			shape[2] = bzg3.size(2);
			auto np = py::array_t<ssize_t,py::array::c_style>(shape);
			size_t nret = bzg3.unsafe_get_map( (ssize_t*)np.request().ptr );
			if (nret != shape[0]*shape[1]*shape[2])
				// I guess nret is smaller, otherwise we probably already encountered a segfault
				throw std::runtime_error("Something has gone horribly wrong with getting the map.");
			return np;
		},
		/*set map*/ [](BrillouinZoneGrid3& bzg3, py::array_t<ssize_t,py::array::c_style> pymap){
			py::buffer_info bi = pymap.request();
			if (bi.ndim != 3) throw std::runtime_error("The mapping must be a 3 dimensional array");
			for (size_t i=0; i<3; ++i) if (bi.shape[i]!=bzg3.size(i))
				throw std::runtime_error("The new map shape must match the old map"); // or we could resize it, but that is more difficult
			if (bzg3.maximum_mapping( (ssize_t*)bi.ptr ) > bzg3.num_data() )
				throw std::runtime_error("The largest integer in the new mapping exceeds the number of data elements.");
			bzg3.unsafe_set_map( (ssize_t*)bi.ptr ); //no error, so this works.
	});
	bzg.def("interpolate_at",[](BrillouinZoneGrid3& bzg3, py::array_t<double,py::array::c_style> pyX, const bool& moveinto){
		py::buffer_info bi = pyX.request();
		if ( bi.shape[bi.ndim-1] !=3 )
			throw std::runtime_error("Interpolation requires one or more 3-vectors");
		ssize_t npts = 1;
		if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
		BrillouinZone b = bzg3.get_brillouinzone();
		Reciprocal lat = b.get_lattice();
		LQVec<double> qv(lat,npts, (double*)bi.ptr ); //memcopy
		if (moveinto){
			LQVec<double> Qv(qv); // second memcopy
			LQVec<int>  tauv(lat,npts); // filled by moveinto
			bool success = b.moveinto(&Qv,&qv,&tauv);
			if (!success)
				throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
		}
		// do the interpolation for each point in qv
		ArrayVector<double> lires = bzg3.linear_interpolate_at(qv);
		// and then make sure we return an numpy array of appropriate size:
		std::vector<ssize_t> outshape;
		for (ssize_t i=0; i < bi.ndim-1; ++i) outshape.push_back(bi.shape[i]);
		if (bzg3.data_ndim() > 1){
			ArrayVector<size_t> data_shape = bzg3.data_shape();
			// the shape of each element is data_shape[1,...,data_ndim-1]
			for (ssize_t i=1; i<data_shape.size(); ++i) outshape.push_back( data_shape.getvalue(i) );
		}
		size_t total_elements = 1;
		for (ssize_t osi : outshape) total_elements *= osi;
		if (total_elements != lires.numel()*lires.size() ){
			std::cout << "outshape: (";
			for (ssize_t osi : outshape) std::cout << osi << "," ;
			std::cout << "\b)" << std::endl;
			printf("Why do we expect %u total elements but only get back %u? (%u × %u)\n",total_elements,lires.numel()*lires.size(),lires.numel(),lires.size());
			throw std::runtime_error("error determining output size");
		}
		auto liout = py::array_t<double,py::array::c_style>(outshape);
		double *rptr = (double*) liout.request().ptr;
		for (size_t i =0; i< lires.size(); i++)
			for (size_t j=0; j< lires.numel(); j++)
				rptr[i*lires.numel()+j] = lires.getvalue(i,j);
		return liout;
	},py::arg("Q"),py::arg("moveinto")=true);

}
