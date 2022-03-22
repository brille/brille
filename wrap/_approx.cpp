#include <pybind11/pybind11.h>
#include "approx_float.hpp"
#include "approx_config.hpp"



void wrap_approx(pybind11::module &m){
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace brille::approx_float;
  m.def("real_space_tolerance", [&](){return config.direct<double>();});
  m.def("reciprocal_space_tolerance", [&](){return config.reciprocal<double>();});
  m.def("real_space_tolerance", [&](double tol){config.direct(tol);});
  m.def("reciprocal_space_tolerance", [&](double tol){config.reciprocal(tol);});

  py::class_<Config> ApproxConfig(m, "ApproxConfig");
  ApproxConfig.def(py::init([](int dig, double dir, double rec){return Config(dig, dir, rec);}), "digits"_a=1000, "real_space_tolerance"_a=1e-10, "reciprocal_space_tolerance"_a=1e-10);
  ApproxConfig.def_property("digits",
    [](const Config & c){return c.digit();},
    [](Config & c, int digits){c.digit(digits); return c.digit();}
  );
  ApproxConfig.def_property("real_space_tolerance",
      [](const Config & c){return c.direct<double>();},
      [](Config & c, double tolerance){c.direct(tolerance); return c.direct<double>();}
  );
  ApproxConfig.def_property("reciprocal_space_tolerance",
      [](const Config & c){return c.reciprocal<double>();},
      [](Config & c, double tolerance){c.reciprocal(tolerance); return c.reciprocal<double>();}
  );
}
