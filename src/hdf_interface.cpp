#include "hdf_interface.hpp"

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <process.h>
std::string brille::hdf_pid_filename(std::string base) {
 return base + std::to_string(_getpid()) + ".h5";
}
#else
#include <unistd.h>
std::string brille::hdf_pid_filename(std::string base) {
  return base + std::to_string(getpid()) + ".h5";
}
#endif



// explicitly define the HighFive create_datatype specialisations rather than
// using the provided macro HIGHFIVE_REGISTER_TYPE
namespace HighFive {
template<> DataType create_datatype<brille::RotatesLike>(){
  return EnumType<brille::RotatesLike>({
      {"Real", brille::RotatesLike::Real},
      {"Reciprocal", brille::RotatesLike::Reciprocal},
      {"Axial", brille::RotatesLike::Axial},
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
      {"inverse_angstrom", brille::LengthUnit::inverse_angstrom}
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