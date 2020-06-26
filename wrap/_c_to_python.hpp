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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>

#include "arrayvector.hpp"
#include "utilities.hpp"

#ifndef __C_TO_PYTHON_H
#define __C_TO_PYTHON_H

namespace py = pybind11;
template<typename T, size_t N> py::array_t<T> sa2np(const std::vector<ssize_t>& sz, const std::array<T,N>& sv){
  // size_t numel = 1;
  // for (ssize_t i: sz) numel *= i;
  size_t numel = signed_to_unsigned<size_t>(std::accumulate(sz.begin(), sz.end(), 1, std::multiplies<ssize_t>()));
  if (N != numel){
    std::string msg = "Inconsistent required shape ( ";
    for (ssize_t i: sz) msg += std::to_string(i) + " ";
    msg += ") and array size " + std::to_string(N);
    throw std::runtime_error(msg);
  }
  auto np = py::array_t<T,py::array::c_style>(sz);
  T *ptr = (T*) np.request().ptr;
  for (size_t i=0; i<numel; ++i) ptr[i] = sv[i];
  return np;
}
template<typename T> py::array_t<T> sv2np(const std::vector<ssize_t>& sz, const std::vector<T>& sv){
  // size_t numel = 1;
  // for (ssize_t i: sz) numel *= i;
  size_t numel = signed_to_unsigned<size_t>(std::accumulate(sz.begin(), sz.end(), 1, std::multiplies<ssize_t>()));
  if (sv.size() != numel){
    std::string msg = "Inconsistent required shape ( ";
    for (ssize_t i: sz) msg += std::to_string(i) + " ";
    msg += ") and vector size " + std::to_string(sv.size());
    throw std::runtime_error(msg);
  }
  auto np = py::array_t<T,py::array::c_style>(sz);
  T *ptr = (T*) np.request().ptr;
  for (size_t i=0; i<numel; ++i) ptr[i] = sv[i];
  return np;
}
template<typename T, size_t N>
py::array_t<T> sva2np(const std::vector<ssize_t>&sz,
                      const std::vector<std::array<T,N>>& sva)
{
  // size_t numel = 1;
  // for (ssize_t i: sz) numel *= i;
  size_t numel = signed_to_unsigned<size_t>(std::accumulate(sz.begin(), sz.end(), 1, std::multiplies<ssize_t>()));
  if (sva.size()*N != numel){
    std::string msg = "Inconsistent required shape ( ";
    for (ssize_t i: sz) msg += std::to_string(i) + " ";
    msg += ") and vector<array<" + std::to_string(N) + "> size ";
    msg += std::to_string(sva.size()) + ".";
    throw std::runtime_error(msg);
  }
  auto np = py::array_t<T, py::array::c_style>(sz);
  T *ptr = (T*) np.request().ptr;
  size_t i=0;
  for (std::array<T,N> array: sva) for (T value: array) ptr[i++]=value;
  return np;
}

template<typename T> py::array_t<T> av2np(const ArrayVector<T>& av){
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
template<typename T> py::array_t<T> av2np_squeeze(const ArrayVector<T>& av){
  std::vector<ssize_t> shape(1); // assume we'll be able to squeeze.
  if (av.size()==1u)
    shape[0] = av.numel();
  else{
    if (av.numel()==1u)
      shape[0] = av.size();
    else
      return av2np(av); // no squeezing possible
  }
  auto np = py::array_t<T,py::array::c_style>(shape);
  T *rptr = (T*) np.request().ptr;
  if (av.size()==1u)
    for (size_t j=0; j<av.numel(); ++j)
      rptr[j] = av.getvalue(0,j);
  else
    for (size_t i=0; i<av.size(); ++i)
      rptr[i] = av.getvalue(i,0);
  return np;
}
/*! \brief Convert an ArrayVector to a numpy.ndarray, using a defined shape

The ArrayVector type is always two-dimensional, but is used by MapGrid3 and
MapGrid4 to represent higher-dimensional data stored at each mapped grid point.
When passing back gridded information to Python, we can re-create the shape
information in the created numpy.ndarray.
@param av The (numel(),size()) ArrayVector
@param inshape A (1,1+array_ndim) ArrayVector consisting of
              (n_arrays, array_dim₀, array_dim₁, …, array_dimₙ)
@param squeeze A flag to indicate if, for av.size()==1, whether to return a
               (1, array_dim₀, array_dim₁, ..., array_dimₙ) or a
               (array_dim₀, array_dim₁, ..., array_dimₙ) numpy.ndarray.
               Or to remove all singleton dimensions from the output.
@returns A pybind wrapped numpy.ndarray with shape
        (av.size(), array_dim₀, array_dim₁, …, array_dimₙ)
*/
template<typename T>
py::array_t<T> av2np_shape(const ArrayVector<T>& av,
                           const ArrayVector<size_t>& inshape,
                           const bool squeeze = false){
  std::vector<ssize_t> outshape;
  if (!(squeeze && av.size()==1))
    outshape.push_back( (ssize_t)av.size() );
  for (size_t i=1; i<inshape.size(); ++i)
    if (!(squeeze && inshape.getvalue(i)==1))
      outshape.push_back( (ssize_t)inshape.getvalue(i) );
  // size_t numel = 1;
  // for (ssize_t osi : outshape) numel *= (size_t)osi;
  size_t numel = signed_to_unsigned<size_t>(std::accumulate(outshape.begin(), outshape.end(), 1, std::multiplies<ssize_t>()));
  if (numel != av.size()*av.numel()){
    return squeeze ? av2np_squeeze(av) : av2np(av);
    // std::string msg = "Expected " + std::to_string(numel)
    //                 + " but only have " + std::to_string(av.numel()*av.size())
    //                 + " (" + std::to_string(av.numel())
    //                 + " × " + std::to_string(av.size()) + ")";
    // throw std::runtime_error(msg);
  }
  auto out = py::array_t<T,py::array::c_style>(outshape);
  T *ptr = (T*)out.request().ptr;
  for (size_t i=0; i<av.size(); i++)
    for (size_t j=0; j<av.numel(); j++)
      ptr[i*av.numel()+j] = av.getvalue(i,j);
  return out;
}
template<typename T>
py::array_t<T> av2np_shape(const ArrayVector<T>& av,
                           const std::vector<size_t>& inshape,
                           const bool squeeze = false){
  std::vector<ssize_t> outshape;
  if (!(squeeze && av.size()==1))
    outshape.push_back( (ssize_t)av.size() );
  for (size_t i=1; i<inshape.size(); ++i)
    if (!(squeeze && inshape[i]==1))
      outshape.push_back( unsigned_to_signed<ssize_t>(inshape[i]) );
  size_t numel = signed_to_unsigned<size_t>(std::accumulate(outshape.begin(), outshape.end(), 1, std::multiplies<ssize_t>()));
  if (numel != av.size()*av.numel()){
    return squeeze ? av2np_squeeze(av) : av2np(av);
  }
  auto out = py::array_t<T,py::array::c_style>(outshape);
  T *ptr = (T*)out.request().ptr;
  for (size_t i=0; i<av.size(); i++)
    for (size_t j=0; j<av.numel(); j++)
      ptr[i*av.numel()+j] = av.getvalue(i,j);
  return out;
}

template<typename T>
std::vector<T> np2vec(py::array_t<T> pyV){
  py::buffer_info vinfo = pyV.request();
  if (vinfo.ndim != 1u)
    throw std::runtime_error("np2vec expects a 1-D input buffer object");
  size_t span = static_cast<size_t>(vinfo.strides[0])/sizeof(T);
  std::vector<T> v;
  size_t vlen = static_cast<size_t>(vinfo.shape[0]);
  v.reserve(vlen);
  T * vptr = (T*) vinfo.ptr;
  for (size_t i=0; i<vlen; ++i) v.push_back(vptr[i*span]);
  return v;  
}

#endif
