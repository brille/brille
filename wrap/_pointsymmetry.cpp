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
#include "_c_to_python.hpp"
#include "pointsymmetry.hpp"
#include "spg_database.hpp"

void wrap_pointsymmetry(pybind11::module & m){
  using pybind11::ssize_t;
  using namespace pybind11::literals;
  using namespace brille;
  pybind11::class_<PointSymmetry> cls(m, "PointSymmetry", R"pbdoc(
  Holds the :math:`3 \times 3` rotation matrices :math:`R` which comprise
  a point group symmetry.

  A point group describes the local symmetry of a lattice point. It contains
  all of the generalised rotations of a :py:class:`~brille._brille.Symmetry`
  with none of its translations.
  )pbdoc");

  // Deprecated: as part of removing the 'Hall number' this initializer will go.
  cls.def(pybind11::init([](int hall, int time_reversal){
    if (PyErr_WarnEx(PyExc_DeprecationWarning,
        "PointSymmetry(Hall_number, time_reversal) is deprecated and will be removed; "
        "use PointSymmetry(Symmetry(...)) or Lattice(...).pointgroup instead. "
        "Time-reversal symmetry is set on BrillouinZone.", 1) < 0)
      throw pybind11::error_already_set();
    return Spacegroup(hall).get_pointgroup_symmetry(time_reversal);}),
    "Hall_number"_a,"time_reversal"_a=0, R"pbdoc(
    Deprecated: use ``PointSymmetry(Symmetry(...))`` or :py:attr:`Lattice.pointgroup`.
    )pbdoc");

  cls.def(pybind11::init([](const Symmetry & sym){
       return PointSymmetry(get_unique_rotations(sym.getallr(), 0));
     }),"Symmetry"_a);

  cls.def_property_readonly("size",&PointSymmetry::size);

  cls.def_property_readonly("W",[](PointSymmetry& ps){
    std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3u, 3u};
    return sva2np(sz, ps.getall());
  });

  cls.def_property_readonly("order",[](PointSymmetry& ps){
    std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size())};
    return sv2np(sz, ps.orders());
  });

  cls.def_property_readonly("isometry",[](PointSymmetry& ps){
    std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size())};
    return sv2np(sz, ps.isometries());
  });

  cls.def_property_readonly("axis",[](PointSymmetry& ps){
    std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3u};
    return sva2np(sz, ps.axes());
  });

  cls.def_property_readonly("generate",&PointSymmetry::generate);

  cls.def_property_readonly("generators",&PointSymmetry::generators);

  cls.def("nfolds",&PointSymmetry::nfolds);
}
