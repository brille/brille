#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>

// the number of nodes we might hold determines what type we need to store
// their indices:
//    bytes   maximum number    type
//      1                255    unsigned short
//      2             65,535    unsigned int
//      4      4,294,967,295    unsigned long
//      8            18×10¹⁸    unsigned long long (aka size_t)
// 65k is not enough. 4B *should* always be sufficient -- each node would
// occupy a fractional volume of ~2×10⁻¹⁰ of the polyhedron, which is overkill.
namespace brille {
  using ind_t = unsigned;
  using shape_t = std::vector<ind_t>;
}

#endif
