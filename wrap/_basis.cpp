#include <pybind11/pybind11.h>

#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <thread>
#include <utility>

#include "_array.hpp"
#include "_c_to_python.hpp"

#include "basis.hpp"

namespace py = pybind11;

void wrap_basis(py::module & m){
  using namespace pybind11::literals;
  using namespace brille;


  py::class_<Basis> cls(m, "Basis");

  cls.def(py::init([](const py::array_t<double>& positions){
    auto pos = np2sva<double,3>(positions);
    return Basis(pos);
  }));
  cls.def(py::init([](const py::array_t<double>& positions, const std::vector<ind_t>& types){
    auto pos = np2sva<double,3>(positions);
    return Basis(pos, types);
  }));

  cls.def_property_readonly("size", &Basis::size);
  cls.def_property_readonly("positions", [](const Basis& b){
    std::vector<size_t> sz{b.size(), 3u};
    return sva2np(sz, b.positions());
  });
  cls.def_property_readonly("types", &Basis::types);
  cls.def("__repr__", &Basis::to_string);
}