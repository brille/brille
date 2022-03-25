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
#include "_nest.hpp"

void wrap_nest(pybind11::module & m){
  using D = double;
  using C = std::complex<double>;
  declare_bznestq<D,D,D>(m,"dd");
  declare_bznestq<D,C,D>(m,"dc");
  declare_bznestq<C,C,D>(m,"cc");
}
