#include <pybind11/pybind11.h>
#include "approx_float.hpp"
#include "approx_config.hpp"



void wrap_approx(pybind11::module &m){
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace brille::approx_float;
  m.def("real_space_tolerance", [&](){return config.direct<double>();},
        R"pbdoc(Return the module-global real space floating point tolerance in angstrom.)pbdoc");
  m.def("reciprocal_space_tolerance", [&](){return config.reciprocal<double>();},
        R"pbdoc(Return the module-global reciprocal space floating point tolerance in inverse angstrom.)pbdoc");
  m.def("real_space_tolerance", [&](double tol){config.direct(tol);}, "tolerance"_a,
        R"pbdoc(Set the module-global real space floating point tolerance in angstrom)pbdoc");
  m.def("reciprocal_space_tolerance", [&](double tol){config.reciprocal(tol);}, "tolerance"_a,
        R"pbdoc(Set the module-global reciprocal space floating point tolerance in angstrom)pbdoc");

  py::class_<Config> ApproxConfig(m, "ApproxConfig",
                                  R"pbdoc(A local set of approximate floating point comparison values)pbdoc");
  ApproxConfig.def(py::init([](int dig, double dir, double rec){return Config(dig, dir, rec);}), "digits"_a=1000, "real_space_tolerance"_a=1e-10, "reciprocal_space_tolerance"_a=1e-10);
  ApproxConfig.def_property("digits",
    [](const Config & c){return c.digit();},
    [](Config & c, int digits){c.digit(digits); return c.digit();},
    R"pbdoc(A multiplier on the machine epsilon, below which two floating point numbers are approximately the same.)pbdoc"
  );
  ApproxConfig.def_property("real_space_tolerance",
      [](const Config & c){return c.direct<double>();},
      [](Config & c, double tolerance){c.direct(tolerance); return c.direct<double>();},
      R"pbdoc(An absolute real space floating point tolerance in angstrom)pbdoc"
  );
  ApproxConfig.def_property("reciprocal_space_tolerance",
      [](const Config & c){return c.reciprocal<double>();},
      [](Config & c, double tolerance){c.reciprocal(tolerance); return c.reciprocal<double>();},
      R"pbdoc(An absolute reciprocal space floating point tolerance in inverse angstrom)pbdoc"
  );
}
