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
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>
#include <tuple>
#include <array>

#include "_interpolator.hpp"

namespace py = pybind11;
namespace br = brille;

void wrap_interpolator(py::module &m){
  using namespace pybind11::literals;
  using namespace brille;
  py::enum_<RotatesLike> enm(m,"RotatesLike",
  R"pbdoc(
    Enumeration indicating how vector and matrix values transform
  )pbdoc"
);
  enm.value("vector", RotatesLike::vector, R"pbdoc(Rotates like a vector)pbdoc");
  enm.value("pseudovector", RotatesLike::pseudovector, R"pbdoc(Rotates like a pseudovector)pbdoc");
  enm.value("Gamma", RotatesLike::Gamma, R"pbdoc(Rotates like a (real space) phonon eigenvector)pbdoc");
}

std::tuple<br::RotatesLike, br::LengthUnit, int, int, std::array<double,3>>
set_check(
  py::array_t<int> pyflg, py::array_t<double> pywght
){
  using namespace brille;
  py::buffer_info bi;
  bi = pyflg.request();
  if (bi.ndim != 1) throw std::runtime_error("flags must be a 1-D array");
  int *intel = (int*) bi.ptr;
  // convert the input integer to a RotatesLike
  RotatesLike rl{RotatesLike::vector};
  if (bi.shape[0] > 0) switch(intel[0]){
    case 2: rl = RotatesLike::Gamma; break;
    case 1: rl = RotatesLike::pseudovector; break;
    case 0: rl = RotatesLike::vector; break;
    default: throw std::runtime_error("Unknown RotatesLike value "+std::to_string(intel[0]));
  }
  // convert input integer to LengthUnit
  LengthUnit lu{LengthUnit::real_lattice};
  if (bi.shape[0] > 1) switch(intel[1]){
    case 4: lu = LengthUnit::reciprocal_lattice; break;
    case 3: lu = LengthUnit::real_lattice; break;
    case 2: lu = LengthUnit::inverse_angstrom; break;
    case 1: lu = LengthUnit::angstrom; break;
    case 0: lu = LengthUnit::none; break;
    default: throw std::runtime_error("Unknown LengthUnit value "+std::to_string(intel[1]));
  }
  // get the cost-function type(s)
  int csf{0}, cvf{0};
  if (bi.shape[0] > 2) csf = intel[2];
  if (bi.shape[0] > 3) cvf = intel[3];
  // copy-over the weight specification
  std::array<double,3> wght{{1,1,1}};
  bi = pywght.request();
  if (bi.ndim != 1) throw std::runtime_error("weights must be a 1-D array");
  auto *dblwght = (double*) bi.ptr;
  for (pybind11::ssize_t i=0; i<bi.shape[0] && i<3; ++i) wght[i] = dblwght[i];
  // tie everything up
  return std::make_tuple(rl, lu, csf, cvf, wght);
}
