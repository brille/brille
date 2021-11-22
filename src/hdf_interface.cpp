#include "hdf_interface.hpp"
#include "rotates.hpp"
#include "trellis_node.hpp"
#include "motion.hpp"

using namespace brille;

HighFive::EnumType<RotatesLike> create_enum_rotateslike(){
  return {
    {"Real", RotatesLike::Real},
    {"Reciprocal", RotatesLike::Reciprocal},
    {"Axial", RotatesLike::Axial},
    {"Gamma", RotatesLike::Gamma}
  };
}

// this macro defines a new function specialization of create_datatype<brille::RotatesLike> in the HIGHFIVE namespace
HIGHFIVE_REGISTER_TYPE(RotatesLike, create_enum_rotateslike)


HighFive::EnumType<NodeType> create_enum_nodetype(){
  return {
    {"null", NodeType::null},
    {"cube", NodeType::cube},
    {"poly", NodeType::poly}
  };
}
HIGHFIVE_REGISTER_TYPE(NodeType, create_enum_nodetype)

HighFive::EnumType<Bravais> create_enum_bravais(){
    // enum class Bravais {_, P, A, B, C, I, F, R};
    return {
            {"_", Bravais::_}, {"P", Bravais::P},
            {"A", Bravais::A}, {"B", Bravais::B}, {"C", Bravais::C},
            {"I", Bravais::I}, {"F", Bravais::F}, {"R", Bravais::R}
    };
}
HIGHFIVE_REGISTER_TYPE(Bravais, create_enum_bravais)

HighFive::EnumType<LengthUnit> create_enum_lengthunit(){
    return {
        {"none", LengthUnit::none},
        {"angstrom", LengthUnit::angstrom},
        {"inverse_angstrom", LengthUnit::inverse_angstrom}
    };
}
HIGHFIVE_REGISTER_TYPE(LengthUnit, create_enum_lengthunit)

template<class T, class R>
HighFive::CompoundType create_compound_Motion(){
    return {
            {"xx", HighFive::AtomicType<T>{}},
            {"xy", HighFive::AtomicType<T>{}},
            {"xz", HighFive::AtomicType<T>{}},
            {"yx", HighFive::AtomicType<T>{}},
            {"yy", HighFive::AtomicType<T>{}},
            {"yz", HighFive::AtomicType<T>{}},
            {"zx", HighFive::AtomicType<T>{}},
            {"zy", HighFive::AtomicType<T>{}},
            {"zz", HighFive::AtomicType<T>{}},
            {"x", HighFive::AtomicType<R>{}},
            {"y", HighFive::AtomicType<R>{}},
            {"z", HighFive::AtomicType<R>{}},
    };
}
//HighFive::CompoundType create_compound_Motion(){
//    return {
//            {"xx", HighFive::AtomicType<int>{}},
//            {"xy", HighFive::AtomicType<int>{}},
//            {"xz", HighFive::AtomicType<int>{}},
//            {"yx", HighFive::AtomicType<int>{}},
//            {"yy", HighFive::AtomicType<int>{}},
//            {"yz", HighFive::AtomicType<int>{}},
//            {"zx", HighFive::AtomicType<int>{}},
//            {"zy", HighFive::AtomicType<int>{}},
//            {"zz", HighFive::AtomicType<int>{}},
//            {"x", HighFive::AtomicType<double>{}},
//            {"y", HighFive::AtomicType<double>{}},
//            {"z", HighFive::AtomicType<double>{}},
//    };
//}
typedef HF_Motion<int,double> HF_Motion_int_double;
HIGHFIVE_REGISTER_TYPE(HF_Motion_int_double, (create_compound_Motion<int,double>))