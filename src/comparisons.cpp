#include <string>
#include "comparisons.hpp"

namespace brille{

template<>
std::string to_string(const cmp& c){
  std::string repr;
  switch(c)
  {
  case cmp::lt:    repr = "<";     break;
  case cmp::gt:    repr = ">";     break;
  case cmp::le:    repr = "<=";    break;
  case cmp::ge:    repr = ">=";    break;
  case cmp::eq:    repr = "==";    break;
  case cmp::nle:   repr = "!<=";   break;
  case cmp::nge:   repr = "!>=";   break;
  case cmp::neq:   repr = "!=";    break;
  case cmp::le_ge: repr = "<=|>="; break;
  default: repr = "unknown comparator";
  }
  return repr;
}

template<>
std::string to_string(const brille::ops& o){
  std::string repr;
  switch(o)
  {
  case ops::plus:  repr = "+";  break;
  case ops::minus: repr = "-";  break;
  case ops::times: repr = "*";  break;
  case ops::rdiv:  repr = "/";  break;
  case ops::ldiv:  repr = "\\"; break;
  default: repr = "unknown operator";
  }
  return repr;
}

}
