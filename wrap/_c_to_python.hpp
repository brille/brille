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
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>

#include "array.hpp"
#include "utilities.hpp"

#ifndef WRAP_BRILLE_C_TO_PYTHON_HPP_
#define WRAP_BRILLE_C_TO_PYTHON_HPP_

namespace py = pybind11;
template<typename T, size_t N> py::array_t<T> sa2np(const std::vector<ssize_t>& sz, const std::array<T,N>& sv){
  // size_t numel = 1;
  // for (ssize_t i: sz) numel *= i;
  size_t numel = brille::utils::s2u<size_t>(std::accumulate(sz.begin(), sz.end(), 1, std::multiplies<ssize_t>()));
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
  size_t numel = brille::utils::s2u<size_t>(std::accumulate(sz.begin(), sz.end(), 1, std::multiplies<ssize_t>()));
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
  size_t numel = brille::utils::s2u<size_t>(std::accumulate(sz.begin(), sz.end(), 1, std::multiplies<ssize_t>()));
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

template<typename T, size_t N>
std::vector<std::array<T,N>>
np2sva(py::array_t<T> np){
  py::buffer_info info = np.request();
  if (info.ndim < 2u)
    throw std::runtime_error("np2sva expects a >1-D input buffer object");
  size_t array_size{1u};
  for (ssize_t i=1; i < info.ndim; ++i) array_size *= info.shape[i];
  if (N != array_size){
    std::string msg = "wrong number of elements beyond first dimension for vector of ";
    msg += std::to_string(N) + " arrays";
    throw std::runtime_error(msg);
  }
  std::vector<std::array<T,N>> v;
  size_t vlen = static_cast<size_t>(info.shape[0]);
  v.reserve(vlen);
  T * ptr = (T*) info.ptr;

  // find the row-ordered spans
  std::vector<size_t> spans(info.ndim, 1u);
  for (ssize_t i=info.ndim-1u; i-->0; ) // first i in loop is ndim-2
    spans[i] = spans[i+1]*info.shape[i+1];
  // check if input np *is* row ordered
  bool roword{true};
  for (ssize_t i=0; i<info.ndim; ++i) roword &= info.strides[i]/sizeof(T) == spans[i];
  if (roword){
    // if it is the copying is straightforward
    for (size_t i=0; i<vlen; ++i){
      std::array<T,N> a;
      for (size_t j=0; j<N; ++j) a[j] = ptr[i*spans[0]+j];
      v.push_back(a);
    }
  } else {
    // we need to worry about subscripted indexing
    for (size_t i=0; i<vlen; ++i){
      std::array<T,N> a;
      // loop over the linear indices within this array (2nd+ dimension[s])
      for (size_t j=0; j<N; ++j){
        // calculate the subscripted indices using the new row-ordered span
        // and the linear index into ptr using the provided strides
        size_t tmp = j; // linear index jᵗʰ array entry
        size_t lin=0; // will be linear index to entry [i,j] of numpy array
        for (ssize_t k=1; k<info.ndim; ++k){
          size_t idx = tmp/spans[k]; // the subscripted index along dimension k
          tmp -= idx*spans[k]; // adjust the linear index
          lin += idx*info.strides[k]/sizeof(T);
        }
        // we skipped the first dimension in the calculation of lin:
        a[j] = ptr[i*info.strides[0]/sizeof(T) + lin];
      }
      v.push_back(a);
    }
  }
  return v;
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
