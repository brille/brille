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
#ifndef BRILLE_BRAVAIS_HPP_
#define BRILLE_BRAVAIS_HPP_
#include <string>
namespace brille {
/*! \brief A Bravais letter indicating a centering of a lattice whose conventional cell is centred.

When the unit cell does not reflec the symmetry of the lattice, it is usual to
refer to a 'conventional' crystallographic basis, (aₛ bₛ cₛ), instead of a
primitive basis, (aₚ bₚ cₚ).
Such a conventional basis has "extra" lattice points added at the centre of the
unit cell, the centre of a face, or the centre of three faces.
The "extra" nodes in the conventional basis are displaced from the origin of the
unit cell by 'centring vectors'. As with any space-spanning basis, any
whole-number linear combination of the conventional basis vectors is a lattice
point but in addition there exist linear combinations xaₛ+ybₛ+zcₛ with at least
two fractional coefficients (x,y,z) that are lattice points as well.

Each conventional basis is ascribed a Bravais letter, which forms part of the
Hermann-Mauguin symbol of a space group.
A subset of the 10 possible Bravais letters is used herein:

| Bravais letter | Centring | Centring vectors |
| --- | --- | --- |
| P | primitve | 0 |
| A | A-face centred | ½bₛ+½cₛ |
| B | B-face centred | ½cₛ+½aₛ |
| C | C-face centred | ½aₛ+½bₛ |
| I | body centred (*Innenzentriert*) | ½aₛ+½bₛ+½cₛ |
| F | all-face centred | ½bₛ+½cₛ, ½cₛ+½aₛ, ½aₛ+½bₛ |
| R | rhombohedrally centred (hexagonal axes) | ⅔aₛ+⅓bₛ+⅓cₛ, ⅓aₛ+⅓bₛ+⅔c |

For further details, see http://reference.iucr.org/dictionary/Centred_lattice
*/
enum class Bravais {_, P, A, B, C, I, F, R};

std::string bravais_string(const Bravais b);
char bravais_letter(const Bravais b);

bool bravais_is_known(const Bravais b);
}
#endif
