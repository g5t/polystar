#include "hdf_interface.hpp"

// explicitly define the HighFive create_datatype specialisations rather than
// using the provided macro HIGHFIVE_REGISTER_TYPE
namespace HighFive {
//template<> DataType create_datatype<polystar::RotatesLike>(){
//  return EnumType<polystar::RotatesLike>({
//      {"vector", polystar::RotatesLike::vector},
//      {"pseudovector", polystar::RotatesLike::pseudovector},
//      {"Gamma", polystar::RotatesLike::Gamma}
//  });
//}
//template<> DataType create_datatype<polystar::NodeType>(){
//  return EnumType<polystar::NodeType>({
//      {"null", polystar::NodeType::null},
//      {"cube", polystar::NodeType::cube},
//      {"poly", polystar::NodeType::poly}
//  });
//}
//template<> DataType create_datatype<polystar::Bravais>(){
//  return EnumType<polystar::Bravais>({
//      {"_", polystar::Bravais::_},
//      {"P", polystar::Bravais::P},
//      {"A", polystar::Bravais::A},
//      {"B", polystar::Bravais::B},
//      {"C", polystar::Bravais::C},
//      {"I", polystar::Bravais::I},
//      {"F", polystar::Bravais::F},
//      {"R", polystar::Bravais::R}
//  });
//}
//template<> DataType create_datatype<polystar::LengthUnit>(){
//  return EnumType<polystar::LengthUnit>({
//      {"none", polystar::LengthUnit::none},
//      {"angstrom", polystar::LengthUnit::angstrom},
//      {"inverse_angstrom", polystar::LengthUnit::inverse_angstrom},
//      {"real_lattice", polystar::LengthUnit::real_lattice},
//      {"reciprocal_lattice", polystar::LengthUnit::reciprocal_lattice},
//  });
//}
template<> DataType create_datatype<polystar::HF_Matrix<int>>(){
  return create_compound_Matrix<int>();
}
  template<> DataType create_datatype<polystar::HF_Matrix<double>>(){
    return create_compound_Matrix<double>();
  }
template<> DataType create_datatype<polystar::HF_Motion<int,double>>(){
  return create_compound_Motion<int,double>();
}
}