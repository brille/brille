#include "hdf_interface.hpp"
#include "rotates.hpp"
#include "trellis_node.hpp"

using namespace brille;

HighFive::EnumType<RotatesLike> create_enum_rotateslike(){
  return {
    {"Real", RotatesLike::Real},
    {"Reciprocal", RotatesLike::Reciprocal},
    {"Axial", RotatesLike::Axial},
    {"Gamma", RotatesLike::Gamma}
  };
}

HIGHFIVE_REGISTER_TYPE(RotatesLike, create_enum_rotateslike)


HighFive::EnumType<NodeType> create_enum_nodetype(){
  return {
    {"null", NodeType::null},
    {"cube", NodeType::cube},
    {"poly", NodeType::poly}
  };
}
HIGHFIVE_REGISTER_TYPE(NodeType, create_enum_nodetype)
