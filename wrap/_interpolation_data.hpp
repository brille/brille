#ifndef __INTERPOLATION_DATA_HPP_
#define __INTERPOLATION_DATA_HPP_

#include <pybind11/pybind11.h>
#include "interpolation_data.hpp"
#include "phonon.hpp"
#include "utilities.hpp"

template<class T>
std::tuple< ArrayVector<T>, std::vector<size_t>, std::array<element_t,3>, RotatesLike >
fill_check(py::array_t<T> pyarray, py::array_t<int> pyel, const size_t count){
  py::buffer_info bi;
  // copy-over the N-D array information
  bi = pyarray.request();
  ArrayVector<T> data( (T*) bi.ptr, bi.shape, bi.strides);
  if (count != data.size()){
    std::string msg;
    msg = "Provided " + std::to_string(data.size()) + " arrays but ";
    msg += std::to_string(count) + " were expected!";
  }
  // and store its shape
  std::vector<size_t> shape;
  for (ssize_t sh: bi.shape) shape.push_back(signed_to_unsigned<size_t>(sh));
  // copy-over the element specification
  std::array<element_t,3> el{{0,0,0}};
  bi = pyel.request();
  if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
  int *intel = (int*) bi.ptr;
  for (ssize_t i=0; i<bi.shape[0] && i<3; ++i)
    el[i] = static_cast<element_t>(intel[i]);
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
  return std::make_tuple(data, shape, el, rl);
}

template<class T>
std::tuple< ArrayVector<T>, std::vector<size_t>, std::array<element_t,3>, RotatesLike, int, int, std::array<double,3>>
fill_check(py::array_t<T> pyarray, py::array_t<int> pyel, py::array_t<double> pywght, const size_t count){
  py::buffer_info bi;
  // copy-over the N-D array information
  bi = pyarray.request();
  ArrayVector<T> data( (T*) bi.ptr, bi.shape, bi.strides);
  if (count != data.size()){
    std::string msg;
    msg = "Provided " + std::to_string(data.size()) + " arrays but ";
    msg += std::to_string(count) + " were expected!";
  }
  // and store its shape
  std::vector<size_t> shape;
  for (ssize_t sh: bi.shape) shape.push_back(signed_to_unsigned<size_t>(sh));
  // copy-over the element specification
  std::array<element_t,3> el{{0,0,0}};
  bi = pyel.request();
  if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
  int *intel = (int*) bi.ptr;
  for (ssize_t i=0; i<bi.shape[0] && i<3; ++i)
    el[i] = static_cast<element_t>(intel[i]);
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
  if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
  double *dblwght = (double*) bi.ptr;
  for (ssize_t i=0; i<bi.shape[0] && i<3; ++i) wght[i] = dblwght[i];
  // tie everything up
  return std::make_tuple(data, shape, el, rl, csf, cvf, wght);
}

template<class T>
py::array_t<T,py::array::c_style>
iid2np(const ArrayVector<T>& av, const InnerInterpolationData<T>& iid, const std::vector<ssize_t>& preshape){
  std::vector<ssize_t> shape;
  std::copy(preshape.begin(), preshape.end(), std::back_inserter(shape));
  if (iid.shape().size()>1){
    const std::vector<size_t>& ds{iid.shape()};
    // the shape of each element is data_shape[1,…,data_ndim-1]
    for (size_t i=1; i<ds.size(); ++i)
      shape.push_back(unsigned_to_signed<ssize_t>(ds[i]));
  }
  size_t prod = signed_to_unsigned<size_t>(std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<ssize_t>()));
  if (prod != av.size()*av.numel()){
    std::string msg;
    msg = "Expected " + std::to_string(prod) + " but only provided ";
    msg += std::to_string(av.numel()*av.size()) + " values to copy ";
    msg += "[( ";
    for (ssize_t i: shape) msg += std::to_string(i) + " ";
    msg += ") vs ("+std::to_string(av.numel())+" x "+std::to_string(av.size())+")]";
    throw std::runtime_error(msg);
  }
  auto out = py::array_t<T,py::array::c_style>(shape);
  T *ptr = (T*) out.request().ptr;
  for (size_t i=0; i<av.size(); ++i) for (size_t j=0; j<av.numel(); ++j)
    ptr[i*av.numel()+j] = av.getvalue(i,j);
  return out;
}

#endif
