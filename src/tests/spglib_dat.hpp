#pragma once

#include "spg_database.hpp"


/* \brief Construct a brille::Symmetry instance from the SPGlib spacegroup database
 *
 * Uses the same encoded operations as the SPGlib database to construct a brille::Symmetry
 * object. This is used to compare the brille::Symmetry object with the one constructed
 * from the Hall symbol.
 */
brille::Symmetry make_spacegroup_symmetry_object(int hall_number);