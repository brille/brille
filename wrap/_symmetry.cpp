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
#include <pybind11/operators.h>

#include "_c_to_python.hpp"
#include "symmetry.hpp"
#include "spg_database.hpp"

/*!
The symmetry operations of a spacegroup
*/
void wrap_symmetry(pybind11::module & m){
  using pybind11::ssize_t;
  using namespace pybind11::literals;
  using namespace brille;
  pybind11::class_<Symmetry> cls(m, "Symmetry", R"pbdoc(
  One or more symmetry operations of a space group.

  A symmetry operation is the combination of a generalised rotation, :math:`W`,
  and translation, :math:`\mathbf{w}`.
  For any position in space, :math:`\mathbf{x}`, the operation transforms 
  :math:`\mathbf{x}` to another equivalent position

  .. math::
    \mathbf{x}' = W \mathbf{x} + \mathbf{w}

  and can equivalently be expressed as :math:`\mathbf{x}' = \mathscr{M}\mathbf{x}`.

  Crystallographic symmetry operations have an order, :math:`o`,
  for which :math:`\mathscr{M}^o = \mathscr{E}` i.e. :math:`o` repeated
  applications of the operation is equivalent to the identity operator.

  A set of symmetry operations can form a group, :math:`\mathbb{G}`, with the property that
  :math:`\mathscr{M}_k = \mathscr{M}_i \mathscr{M}_j` with :math:`\mathscr{M}_i,\mathscr{M}_j,\mathscr{M}_k \in \mathbb{G}`.

  This class can be used to hold any number of related symmetry operators, and to generate all spacegroup operators from those stored.

  Parameters
  ----------
  hall : int
    The integer Hall number for the desired space group operations [[deprecated]].
  W : arraylike, int
    The generalised rotation (matrix) part of the symmetry operator(s)
  w : arraylike, float
    The translation (vector) part of the symmetry operator(s)
  cifxyz : str
    The symmetry operator(s) encoded in CIF xyz format

  Note
  ----
  The overloaded forms of ``__init__`` take one of **hall**, (**W**, **w**), *or* **cifxyz**.
  )pbdoc");


  cls.def(pybind11::init([](int hall){return Spacegroup(hall).get_spacegroup_symmetry();}),"Hall number"_a);

  cls.def(pybind11::init([](py::array_t<int> pyW, py::array_t<double> pyw){
    std::vector<std::array<int,9>> Ws = np2sva<int,9>(pyW);
    std::vector<std::array<double,3>> ws = np2sva<double,3>(pyw);
    if (Ws.size() != ws.size())
      throw std::runtime_error("Equal number of matrices and vectors required");
    Symmetry::Motions mots;
    mots.reserve(Ws.size());
    for (size_t i=0; i<Ws.size(); ++i)
      mots.push_back(Motion<int,double>(Ws[i], ws[i]));
    Symmetry s(mots);
    return s;
  }),"W"_a,"w"_a);

  cls.def(pybind11::init([](const std::string& str){
    Symmetry s;
    s.from_ascii(str);
    return s;
  }),"CIF xyz string"_a);

  cls.def_property_readonly("size",&Symmetry::size);

  cls.def("__len__", &Symmetry::size);

  cls.def_property_readonly("W",[](Symmetry& ps){
    std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3u, 3u};
    return sva2np(sz, ps.getallr());
  });

  cls.def_property_readonly("w",[](Symmetry& ps){
    std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3u};
    return sva2np(sz, ps.getallt());
  });

  cls.def("generate",&Symmetry::generate);

  cls.def("generators",&Symmetry::generators);

  cls.def_property_readonly("centring",&Symmetry::getcentring);

  cls.def(pybind11::self == pybind11::self);
}
