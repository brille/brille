#ifndef __INTERPOLATION_DATA_HPP_
#define __INTERPOLATION_DATA_HPP_

#include <pybind11/pybind11.h>
#include "interpolation_data.hpp"
#include "phonon.hpp"
#include "utilities.hpp"
#include "_array.hpp"

namespace py = pybind11;

template<class T>
std::tuple< brille::Array<T>, std::array<brille::ind_t,3>, RotatesLike >
fill_check(py::array_t<T> pyarray, py::array_t<int> pyel, const size_t count){
  // this will probably be a problem if I don't figure out a way to add to the python reference count
  brille::Array<T> data = brille::py2a(pyarray);
  if (count != data.size(0)){
    std::string msg;
    msg = "Provided " + std::to_string(data.size(0)) + " arrays but ";
    msg += std::to_string(count) + " were expected!";
  }
  // InnerInterpolationData requires data arrays of atleast 2D:
  if (data.ndim() < 2){
    if (!data.is_contiguous()) data = data.contiguous_copy();
    data.reshape(brille::shape_t{data.size(0),1});
  }
  // copy-over the element specification
  std::array<brille::ind_t,3> el{{0,0,0}};
  py::buffer_info bi = pyel.request();
  if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
  int *intel = (int*) bi.ptr;
  for (ssize_t i=0; i<bi.shape[0] && i<3; ++i)
    el[i] = static_cast<brille::ind_t>(intel[i]);
  // convert the input integer to a RotatesLike
  RotatesLike rl{RotatesLike::Real};
  if (bi.shape[0] > 3) switch(intel[3]){
    case 3: rl = RotatesLike::Gamma; break;
    case 2: rl = RotatesLike::Axial; break;
    case 1: rl = RotatesLike::Reciprocal; break;
    case 0: rl = RotatesLike::Real; break;
    default: throw std::runtime_error("Unknown RotatesLike value "+std::to_string(intel[3]));
  }
  // tie everything up
  return std::make_tuple(data, el, rl);
}

template<class T>
std::tuple< brille::Array<T>, std::array<brille::ind_t,3>, RotatesLike, int, int, std::array<double,3> >
fill_check(py::array_t<T> pyarray, py::array_t<int> pyel, py::array_t<double> pywght, const size_t count){
  // copy-over the N-D array information
  // this will probably be a problem if I don't figure out a way to add to the python reference count
  brille::Array<T> data = brille::py2a(pyarray);
  if (count != data.size(0)){
    std::string msg;
    msg = "Provided " + std::to_string(data.size(0)) + " arrays but ";
    msg += std::to_string(count) + " were expected!";
  }
  // InnerInterpolationData requires data arrays of atleast 2D:
  if (data.ndim() < 2){
    if (!data.is_contiguous()) data = data.contiguous_copy();
    data.reshape(brille::shape_t{data.size(0),1});
  }
  // copy-over the element specification
  std::array<brille::ind_t,3> el{{0,0,0}};
  py::buffer_info bi = pyel.request();
  if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
  int *intel = (int*) bi.ptr;
  for (ssize_t i=0; i<bi.shape[0] && i<3; ++i)
    el[i] = static_cast<brille::ind_t>(intel[i]);
  // convert the input integer to a RotatesLike
  RotatesLike rl{RotatesLike::Real};
  if (bi.shape[0] > 3) switch(intel[3]){
    case 3: rl = RotatesLike::Gamma; break;
    case 2: rl = RotatesLike::Axial; break;
    case 1: rl = RotatesLike::Reciprocal; break;
    case 0: rl = RotatesLike::Real; break;
    default: throw std::runtime_error("Unknown RotatesLike value "+std::to_string(intel[3]));
  }
  // get the cost-function type(s)
  int csf{0}, cvf{0};
  if (bi.shape[0] > 4) csf = intel[4];
  if (bi.shape[0] > 5) cvf = intel[5];
  // copy-over the weight specification
  std::array<double,3> wght{{1,1,1}};
  bi = pywght.request();
  if (bi.ndim != 1) throw std::runtime_error("weights must be a 1-D array");
  double *dblwght = (double*) bi.ptr;
  for (ssize_t i=0; i<bi.shape[0] && i<3; ++i) wght[i] = dblwght[i];
  // tie everything up
  return std::make_tuple(data, el, rl, csf, cvf, wght);
}

std::tuple<RotatesLike, int, int, std::array<double,3>>
set_check(py::array_t<int> pyflg, py::array_t<double> pywght);

#endif
