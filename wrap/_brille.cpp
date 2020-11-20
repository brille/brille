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
#include "version.hpp"

#include <pybind11/pybind11.h>

void wrap_bravais(pybind11::module &);
void wrap_brillouinzone(pybind11::module &);
void wrap_debug(pybind11::module &);
void wrap_hallsymbol(pybind11::module &);
void wrap_lattice(pybind11::module &);
void wrap_mesh(pybind11::module &);
void wrap_nest(pybind11::module &);
void wrap_pointgroup(pybind11::module &);
void wrap_pointsymmetry(pybind11::module &);
void wrap_polyhedron(pybind11::module &);
void wrap_primitivetransform(pybind11::module &);
void wrap_sortingstatus(pybind11::module &);
void wrap_spacegroup(pybind11::module &);
void wrap_symmetry(pybind11::module &);
void wrap_trellis(pybind11::module &);

void wrap_version(pybind11::module & m){
  using namespace brille::version;
  m.attr("__version__") = version_number;
  std::string v = version_number;
  if (!std::string(git_revision).empty()){
    v += "+" + std::string(git_branch);
    v += "." + std::string(git_revision).substr(0,7);
  }
  m.attr("version") = v;
  m.attr("git_revision") = git_revision;
  m.attr("build_datetime") = build_datetime;
  m.attr("build_hostname") = build_hostname;
}

PYBIND11_MODULE(_brille, m){
  m.doc() = R"pbdoc(
    pybind11 module ``brille._brille``
    ---------------------------
    This module provides the interface to the C++ library.

    All of the symbols defined within ``_brille`` are imported by ``brille`` to
    make using them that much easier. If in doubt, the interfaced classes can be
    accessed via their submodule syntax ``brille._brille``.

    .. currentmodule:: brille._brille

    .. autosummary::
      :toctree: _generate

  )pbdoc";
  wrap_version(m);
  wrap_lattice(m);
  wrap_brillouinzone(m);
  wrap_mesh(m);
  wrap_trellis(m);
  wrap_nest(m);
  wrap_primitivetransform(m);
  wrap_spacegroup(m);
  wrap_pointgroup(m);
  wrap_sortingstatus(m);
  wrap_symmetry(m);
  wrap_pointsymmetry(m);
  wrap_polyhedron(m);
  wrap_hallsymbol(m);
  wrap_bravais(m);
  wrap_debug(m);
}
