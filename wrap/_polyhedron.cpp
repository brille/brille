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
#include "_polyhedron.hpp"
void wrap_polyhedron(pybind11::module &m){
  using namespace brille;
  using namespace brille::lattice;
  using namespace brille::polyhedron;
  namespace py = pybind11;

  py::class_<Poly<double,bArray>> bare(m, "Polyhedron");
  define_polyhedron_inits<double>(bare);
  define_polyhedron<double>(bare);

  py::class_<Poly<double,LVec>> lvec(m, "LPolyhedron");
  define_polyhedron<double>(lvec);
  define_polyhedron_lvec<double>(lvec);
}
