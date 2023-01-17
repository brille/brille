#include "vertex_map_set.h"
#include <iostream>

std::ostream & operator<<(std::ostream & os, brille::MapVertexType t){
  using namespace brille;
  switch(t){
  case MapVertexType::Pristine: os << 'P'; break;
  case MapVertexType::Appended: os << 'A'; break;
  case MapVertexType::SecondPristine: os << "sP"; break;
  case MapVertexType::SecondAppended: os << "sA"; break;
  }
  return os;
}