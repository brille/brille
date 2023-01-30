
#include <pybind11/pybind11.h>
#include "interpolatordual.hpp"
#include "phonon.hpp"
#include "utilities.hpp"
#include "_array.hpp"

#ifndef WRAP_BRILLE_INTERPOLATOR_HPP_
#define WRAP_BRILLE_INTERPOLATOR_HPP_
namespace py = pybind11;
namespace br = brille;

template<class T>
br::Interpolator<T>
fill_check(py::array_t<T> pyarray, py::array_t<int> pyel, const size_t count){
  using namespace brille;
  br::Array<T> data = br::py2a(pyarray);
  if (count != data.size(0)){
    std::string msg;
    msg = "Provided " + std::to_string(data.size(0)) + " arrays but ";
    msg += std::to_string(count) + " were expected!";
  }
  // copy-over the element specification
  std::array<br::ind_t,3> el{{0,0,0}};
  py::buffer_info bi = pyel.request();
  if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
  int *intel = (int*) bi.ptr;
  for (pybind11::ssize_t i=0; i<bi.shape[0] && i<3; ++i)
    el[i] = static_cast<br::ind_t>(intel[i]);
  // convert the input integer to a RotatesLike
  RotatesLike rl{RotatesLike::vector};
  if (bi.shape[0] > 3) switch(intel[3]){
    case 2: rl = RotatesLike::Gamma; break;
    case 1: rl = RotatesLike::pseudovector; break;
    case 0: rl = RotatesLike::vector; break;
    default: throw std::runtime_error("Unknown RotatesLike value "+std::to_string(intel[3]));
  }
  // convert input integer to LengthUnit
  LengthUnit lu{LengthUnit::real_lattice};
  if (bi.shape[0] > 4) switch(intel[4]){
    case 4: lu = LengthUnit::reciprocal_lattice; break;
    case 3: lu = LengthUnit::real_lattice; break;
    case 2: lu = LengthUnit::inverse_angstrom; break;
    case 1: lu = LengthUnit::angstrom; break;
    case 0: lu = LengthUnit::none; break;
    default: throw std::runtime_error("Unknown LengthUnit value "+std::to_string(intel[4]));
  }
  // tie everything up
  // return std::make_tuple(data, el, rl);
  return Interpolator(data, el, rl, lu);
}

template<class T>
br::Interpolator<T>
// std::tuple< br::Array<T>, std::array<br::ind_t,3>, RotatesLike, int, int, std::array<double,3> >
fill_check(py::array_t<T> pyarray, py::array_t<int> pyel, py::array_t<double> pywght, const size_t count){
  using namespace brille;
  // copy-over the N-D array information
  // this will probably be a problem if I don't figure out a way to add to the python reference count
  br::Array<T> data = br::py2a(pyarray);
  if (count != data.size(0)){
    std::string msg;
    msg = "Provided " + std::to_string(data.size(0)) + " arrays but ";
    msg += std::to_string(count) + " were expected!";
  }
  // InnerInterpolationData requires data arrays of atleast 2D:
  if (data.ndim() < 2){
    if (!data.is_contiguous()) data = data.contiguous_copy();
    data.reshape(br::shape_t{data.size(0),1});
  }
  // copy-over the element specification
  std::array<br::ind_t,3> el{{0,0,0}};
  py::buffer_info bi = pyel.request();
  if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
  int *intel = (int*) bi.ptr;
  for (pybind11::ssize_t i=0; i<bi.shape[0] && i<3; ++i)
    el[i] = static_cast<br::ind_t>(intel[i]);
  // convert the input integer to a RotatesLike
  RotatesLike rl{RotatesLike::vector};
  if (bi.shape[0] > 3) switch(intel[3]){
    case 2: rl = RotatesLike::Gamma; break;
    case 1: rl = RotatesLike::pseudovector; break;
    case 0: rl = RotatesLike::vector; break;
    default: throw std::runtime_error("Unknown RotatesLike value "+std::to_string(intel[3]));
  }
  // convert input integer to LengthUnit
  LengthUnit lu{LengthUnit::real_lattice};
  if (bi.shape[0] > 4) switch(intel[4]){
    case 4: lu = LengthUnit::reciprocal_lattice; break;
    case 3: lu = LengthUnit::real_lattice; break;
    case 2: lu = LengthUnit::inverse_angstrom; break;
    case 1: lu = LengthUnit::angstrom; break;
    case 0: lu = LengthUnit::none; break;
    default: throw std::runtime_error("Unknown LengthUnit value "+std::to_string(intel[4]));
  }
  // get the cost-function type(s)
  int csf{0}, cvf{0};
  if (bi.shape[0] > 5) csf = intel[5];
  if (bi.shape[0] > 6) cvf = intel[6];
  // copy-over the weight specification
  std::array<double,3> wght{{1,1,1}};
  bi = pywght.request();
  if (bi.ndim != 1) throw std::runtime_error("weights must be a 1-D array");
  auto *dblwght = (double*) bi.ptr;
  for (pybind11::ssize_t i=0; i<bi.shape[0] && i<3; ++i) wght[i] = dblwght[i];
  // tie everything up
  // return std::make_tuple(data, el, rl, csf, cvf, wght);
  return Interpolator(data, el, rl, lu, csf, cvf, wght);
}

std::tuple<br::RotatesLike, br::LengthUnit, int, int, std::array<double,3>>
set_check(py::array_t<int> pyflg, py::array_t<double> pywght);

#endif
