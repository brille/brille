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
#include "bravais.hpp"

using namespace brille;

std::string brille::bravais_string(const Bravais b){
  std::string repr;
  switch (b){
    case Bravais::_: repr = "centering error";        break;
    case Bravais::P: repr = "primitive";              break;
    case Bravais::A: repr = "A-face centred";         break;
    case Bravais::B: repr = "B-face centred";         break;
    case Bravais::C: repr = "C-face centred";         break;
    case Bravais::I: repr = "body centred";           break;
    case Bravais::F: repr = "face centred";           break;
    case Bravais::R: repr = "rhombohedrally centred"; break;
  }
  return repr;
}

char brille::bravais_letter(const Bravais b){
  switch (b){
    case Bravais::_: return '!';
    case Bravais::P: return 'P';
    case Bravais::A: return 'A';
    case Bravais::B: return 'B';
    case Bravais::C: return 'C';
    case Bravais::I: return 'I';
    case Bravais::F: return 'F';
    case Bravais::R: return 'R';
  }
  return '\0';
}

bool brille::bravais_is_known(const Bravais b){
  switch (b){
    case Bravais::P:
    case Bravais::A:
    case Bravais::B:
    case Bravais::C:
    case Bravais::I:
    case Bravais::F:
    case Bravais::R:
      return true;
    default:
      return false;
  }
  return false;
}
