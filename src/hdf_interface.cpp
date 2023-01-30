#include "hdf_interface.hpp"

// explicitly define the HighFive create_datatype specialisations rather than
// using the provided macro HIGHFIVE_REGISTER_TYPE
namespace HighFive {
template<> DataType create_datatype<brille::RotatesLike>(){
  return EnumType<brille::RotatesLike>({
      {"vector", brille::RotatesLike::vector},
      {"pseudovector", brille::RotatesLike::pseudovector},
      {"Gamma", brille::RotatesLike::Gamma}
  });
}
template<> DataType create_datatype<brille::NodeType>(){
  return EnumType<brille::NodeType>({
      {"null", brille::NodeType::null},
      {"cube", brille::NodeType::cube},
      {"poly", brille::NodeType::poly}
  });
}
template<> DataType create_datatype<brille::Bravais>(){
  return EnumType<brille::Bravais>({
      {"_", brille::Bravais::_},
      {"P", brille::Bravais::P},
      {"A", brille::Bravais::A},
      {"B", brille::Bravais::B},
      {"C", brille::Bravais::C},
      {"I", brille::Bravais::I},
      {"F", brille::Bravais::F},
      {"R", brille::Bravais::R}
  });
}
template<> DataType create_datatype<brille::LengthUnit>(){
  return EnumType<brille::LengthUnit>({
      {"none", brille::LengthUnit::none},
      {"angstrom", brille::LengthUnit::angstrom},
      {"inverse_angstrom", brille::LengthUnit::inverse_angstrom},
      {"real_lattice", brille::LengthUnit::real_lattice},
      {"reciprocal_lattice", brille::LengthUnit::reciprocal_lattice},
  });
}
template<> DataType create_datatype<brille::HF_Matrix<int>>(){
  return create_compound_Matrix<int>();
}
  template<> DataType create_datatype<brille::HF_Matrix<double>>(){
    return create_compound_Matrix<double>();
  }
template<> DataType create_datatype<brille::HF_Motion<int,double>>(){
  return create_compound_Motion<int,double>();
}
}