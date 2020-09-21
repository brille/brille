/* Copyright 2020 Greg Tucker
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
#ifndef WRAP_ARRAY_HPP
#define WRAP_ARRAY_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>

#include "array.hpp"

namespace brille {
  template<class T>
  pybind11::array_t<T> a2py(brille::Array<T>& a){
    // share an Array with Python
    // construct a new Array<T> using the same underlying heap data
    std::unique_ptr<brille::Array<T>> aptr = std::make_unique<brille::Array<T>>(brille::Array<T>(a));
    auto capsule = pybind11::capsule(aptr.get(), [](void *p) { std::unique_ptr<brille::Array<T>>(reinterpret_cast<brille::Array<T>*>(p)); });
    aptr.release();
    return pybind11::array_t<T>(a.shape(), a.cstride(), a.data(), capsule);
  }

  // inspired by https://github.com/pybind/pybind11/issues/1042#issuecomment-647147819
  template<class T>
  pybind11::array_t<T> a2py(brille::Array<T>&& a){
    // move an Array to Python.
  	// Ref: https://stackoverflow.com/questions/54876346/pybind11-and-stdvector-how-to-free-data-using-capsules
    auto* aptr = new brille::Array<T>(std::move(a));
    // At this point, transferToHeapGetRawPtr is a raw pointer to an object on the heap.
    // No unique_ptr or shared_ptr, it will have to be freed with delete to avoid a memory leak.
    auto capsule = pybind11::capsule(aptr, [](void *toFree){ delete static_cast<brille::Array<T>*>(toFree); });
    return pybind11::array_t<T>(a.shape(), a.cstride(), aptr->data(), capsule);
  }

  template<class T>
  brille::Array<T> py2a(pybind11::array_t<T> pya){
    pybind11::buffer_info info = pya.request();
    std::vector<ind_t> shape, stride;
    for (ssize_t i=0; i<info.ndim; ++i){
      shape.push_back(static_cast<ind_t>(info.shape[i]));
      stride.push_back(static_cast<ind_t>(info.strides[i]/sizeof(T)));
    }
    T* ptr = (T*) info.ptr;
    ind_t num = info.size;
    return Array<T>(ptr, num, false /*brille doesn't own ptr*/, shape, stride);
  }
}

template<class T>
void declare_array(pybind11::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using Class = brille::Array<T>;
  using ind_t = brille::ind_t;
  std::string pyclass_name = std::string("Array")+typestr;
  pybind11::class_<Class> cls(m, pyclass_name.c_str(), pybind11::buffer_protocol(), pybind11::dynamic_attr());
  //buffer_info
  cls.def_buffer([](Class &cobj) -> pybind11::buffer_info {
    // return brille::a2bi(cobj);
    return pybind11::buffer_info(cobj.data(), cobj.shape(), cobj.cstride(), cobj.ismutable());
  });
  //initializer(s):
  cls.def(pybind11::init<const std::vector<ind_t>&>(),"shape"_a);
  cls.def(pybind11::init<const std::vector<ind_t>&,T>(),"shape"_a, "init"_a);
  cls.def(pybind11::init<const std::vector<ind_t>&,const std::vector<ind_t>&>(), "shape"_a, "stride"_a);
  cls.def(pybind11::init<const std::vector<ind_t>&,const std::vector<ind_t>&,T>(), "shape"_a, "stride"_a, "init"_a);

  cls.def(pybind11::init([](pybind11::array_t<T> buffer_obj){
    return brille::py2a<T>(buffer_obj);
  }),"BufferObject"_a);

  cls.def("to_string",&Class::to_string);
  cls.def("__repr__",&Class::to_string);
}

#endif
