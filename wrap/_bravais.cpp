/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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
#include "bravais.hpp"

namespace py = pybind11;

void wrap_bravais(py::module &m){
  using namespace brille;
  py::enum_<Bravais>(m, "Bravais", R"pbdoc(
  A Bravais letter indicating the centering of a lattice

  When the unit cell does not reflect the symmetry of the lattice, it is usual
  to refer to a 'conventional' crystallographic basis,
  :math:`(\mathbf{a}_s\,\mathbf{b}_s\,\mathbf{c}_s)`, instead of
  a primitive basis, :math:`(\mathbf{a}_p\,\mathbf{b}_p\,\mathbf{c}_p)`.
  Such a conventional basis has 'extra' lattice points added at the centre of
  the unit cell, the centre of a face, or the centre of three faces.
  The 'extra' nodes in the conventional basis are displaced from the origin of
  the unit cell by 'centring vectors'. As with any space-spanning basis, any
  whole-number linear combination of the conventional basis vectors is a lattice
  point but in addition there exist linear combinations
  :math:`x\mathbf{a}_s+y\mathbf{b}_s+z\mathbf{c}_s` with at least
  two fractional coefficients :math:`(x,y,z)` that are lattice points as well.

  Each conventional basis is ascribed a Bravais letter, which forms part of the
  Hermann-Mauguin symbol of a space group.
  A subset of the 10 possible Bravais letters is used herein:

  +----------------+-----------------------------------------+-----------------------------------------------------------+
  | Bravais letter | Centring                                | Centring vectors                                          |
  +================+=========================================+===========================================================+
  |       P        | primitive                               | :math:`\mathbf{0}`                                        |
  +----------------+-----------------------------------------+-----------------------------------------------------------+
  |       A        | A-face centred                          | :math:`\frac{\mathbf{b}_s+\mathbf{c}_s}{2}`               |
  +----------------+-----------------------------------------+-----------------------------------------------------------+
  |       B        | B-face centred                          | :math:`\frac{\mathbf{c}_s+\mathbf{a}_s}{2}`               |
  +----------------+-----------------------------------------+-----------------------------------------------------------+
  |       C        | C-face centred                          | :math:`\frac{\mathbf{a}_s+\mathbf{b}_s}{2}`               |
  +----------------+-----------------------------------------+-----------------------------------------------------------+
  |       I        | body centred (*Innenzentriert*)         | :math:`\frac{\mathbf{a}_s+\mathbf{b}_s+\mathbf{c}_s}{2}`  |
  +----------------+-----------------------------------------+-----------------------------------------------------------+
  |       F        | all-face centred                        | :math:`\frac{\mathbf{b}_s+\mathbf{c}_s}{2}`,              |
  |                |                                         | :math:`\frac{\mathbf{c}_s+\mathbf{a}_s}{2}`,              |
  |                |                                         | :math:`\frac{\mathbf{a}_s+\mathbf{b}_s}{2}`               |
  +----------------+-----------------------------------------+-----------------------------------------------------------+
  |       R        | rhombohedrally centred (hexagonal axes) | :math:`\frac{2\mathbf{a}_s+\mathbf{b}_s+\mathbf{c}_s}{3}` |
  |                |                                         | :math:`\frac{\mathbf{a}_s+2\mathbf{b}_s+2\mathbf{c}_s}{3}`|
  +----------------+-----------------------------------------+-----------------------------------------------------------+

  For further details, see the `IUCr Online Dictionary of Crystallography`__.

  .. _website: http://reference.iucr.org/dictionary/Centred_lattice
  __ website_
  )pbdoc")
    .value("invalid", Bravais::_)
    .value("P", Bravais::P, R"pbdoc(primitive)pbdoc")
    .value("A", Bravais::A, R"pbdoc(A-face centred)pbdoc")
    .value("B", Bravais::B, R"pbdoc(B-face centred)pbdoc")
    .value("C", Bravais::C, R"pbdoc(C-face centred)pbdoc")
    .value("I", Bravais::I, R"pbdoc(body-centred)pbdoc")
    .value("F", Bravais::F, R"pbdoc(face centred)pbdoc")
    .value("R", Bravais::R, R"pbdoc(rhombohedrally centred)pbdoc");
}
