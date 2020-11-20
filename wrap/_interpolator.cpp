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

    Values
    ------
    RotatesLike::Real
      transform like a real space vector or matrix
    RotatesLike::Reciprocal
      transform like a reciprocal space vector or matrix
    RotatesLike::Axis
      transform like an real space axial vector
    RotatesLike::Gamma
      transform like a (real space) phonon eigenvector
  )pbdoc"
);
  enm.value("Real", RotatesLike::Real);
  enm.value("Reciprocal", RotatesLike::Reciprocal);
  enm.value("Axial", RotatesLike::Axial);
  enm.value("Gamma", RotatesLike::Gamma);
}

std::tuple<br::RotatesLike, int, int, std::array<double,3>>
set_check(
  py::array_t<int> pyflg, py::array_t<double> pywght
){
  using namespace brille;
  py::buffer_info bi;
  bi = pyflg.request();
  if (bi.ndim != 1) throw std::runtime_error("flags must be a 1-D array");
  int *intel = (int*) bi.ptr;
  // convert the input integer to a RotatesLike
  RotatesLike rl{RotatesLike::Real};
  if (bi.shape[0] > 0) switch(intel[0]){
    case 3: rl = RotatesLike::Gamma; break;
    case 2: rl = RotatesLike::Axial; break;
    case 1: rl = RotatesLike::Reciprocal; break;
    case 0: rl = RotatesLike::Real; break;
    default: throw std::runtime_error("Unknown RotatesLike value "+std::to_string(intel[3]));
  }
  // get the cost-function type(s)
  int csf{0}, cvf{0};
  if (bi.shape[0] > 1) csf = intel[1];
  if (bi.shape[0] > 2) cvf = intel[2];
  // copy-over the weight specification
  std::array<double,3> wght{{1,1,1}};
  bi = pywght.request();
  if (bi.ndim != 1) throw std::runtime_error("weights must be a 1-D array");
  double *dblwght = (double*) bi.ptr;
  for (ssize_t i=0; i<bi.shape[0] && i<3; ++i) wght[i] = dblwght[i];
  // tie everything up
  return std::make_tuple(rl, csf, cvf, wght);
}
