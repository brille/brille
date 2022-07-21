//
// Created by g on 19/11/2021.
//

#ifndef BRILLE_HDF_INTERFACE_HPP
#define BRILLE_HDF_INTERFACE_HPP
#include "enums.hpp"
#ifdef USE_HIGHFIVE

#include <map>
#include <vector>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>



const std::string LENGTH_ATTR_NAME("length");

namespace brille {

  template<class T, class... V>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, T>, HighFive::Group>
  overwrite_group(T& obj, const std::string& entry, V... args){
      if (obj.exist(entry)) obj.unlink(entry);
      return obj.createGroup(entry, args...);
  }
  template<class T, class... V>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, T>, HighFive::DataSet>
    overwrite_data(T& obj, const std::string& entry, V... args){
      if (obj.exist(entry)) obj.unlink(entry);
      return obj.createDataSet(entry, args...);
  }

  template<class T, class R>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
  lists_to_hdf(const std::vector<std::vector<T>>& lists, R& obj, const std::string& entry){
    using namespace HighFive;
    // if we got this far and are supposed to write into group 'entry',
    // we *must* remove any existing group first:
    if (obj.exist(entry)) obj.unlink(entry);
    // explicitly create the group so that we can add an attribute of its size
    auto group = obj.createGroup(entry);
    group.createAttribute(LENGTH_ATTR_NAME, lists.size());
    size_t i=0;
    for (const auto& list: lists) group.createDataSet(std::to_string(i++), list);
    return true;
  }

  template<class T>
  bool lists_to_hdf(const std::vector<std::vector<T>>& lists, const std::string& filename, const std::string& entry, const unsigned perm=HighFive::File::OpenOrCreate){
    using namespace HighFive;
    File file(filename, perm);
    return lists_to_hdf(lists, file, entry);
  }

  template<class T, class R>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, std::vector<std::vector<T>>>
  lists_from_hdf(R& obj, const std::string& entry){
    auto group = obj.getGroup(entry);
    auto attr = group.getAttribute(LENGTH_ATTR_NAME);
    size_t no;
    attr.read(no);
    std::vector<std::vector<T>> out(no);
    for (size_t i=0; i<no; ++i){
        auto dataset = group.getDataSet(std::to_string(i));
        dataset.read(out[i]);
    }
    return out;
  }

  template<class T>
  std::vector<std::vector<T>> lists_from_hdf(const std::string& filename, const std::string& entry){
      using namespace HighFive;
      File file(filename, File::ReadOnly);
      return lists_from_hdf<T>(file, entry);
  }

  template<class Key, class Value, class R>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
  map_to_hdf(const std::map<Key, Value>& map, R& obj, const std::string& entry){
    using namespace HighFive;
    if (obj.exist(entry)) obj.unlink(entry);
    auto group = obj.createGroup(entry);
    std::vector<Key> keys; keys.reserve(map.size());
    std::vector<Value> values; values.reserve(map.size());
    for (const auto & [key, value]: map){
      keys.push_back(key);
      values.push_back(value);
    }
    group.createDataSet("keys", keys);
    group.createDataSet("values", values);
    return true;
  }
  template<class Key, class Value, class R>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, std::map<Key, Value>>
  map_from_hdf(R& obj, const std::string& entry){
    using namespace HighFive;
    auto group = obj.getGroup(entry);
    std::vector<Key> keys;
    std::vector<Value> values;
    group.getDataSet("keys").read(keys);
    group.getDataSet("values").read(values);
    if (keys.size() != values.size())
      throw std::runtime_error("HDF5 file std::map error");
    std::map<Key, Value> map;
    for (size_t i=0; i<keys.size(); ++i) map.emplace(keys[i], values[i]);
    return map;
  }

    template<class T, class R, size_t N>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
    lists_to_hdf(const std::array<std::vector<T>, N>& lists, R& obj, const std::string& entry){
        using namespace HighFive;
        auto group = overwrite_group(obj, entry);
        group.createAttribute(LENGTH_ATTR_NAME, lists.size());
        size_t i=0;
        for (const auto& list: lists) group.createDataSet(std::to_string(i++), list);
        return true;
    }

  // Special simple version of Motion for HighFive HDF5 interface
  template<class T, class R>
  struct HF_Motion {
    T xx, xy, xz, yx, yy, yz, zx, zy, zz;
    R x, y, z;
  public:
    explicit HF_Motion() = default;
    HF_Motion(const std::array<T, 9>& m, const std::array<R, 3>& v)
    : xx(m[0]), xy(m[1]), xz(m[2]), yx(m[3]), yy(m[4]), yz(m[5]), zx(m[6]), zy(m[7]), zz(m[8]),
      x(v[0]), y(v[1]), z(v[2]) {}
    [[nodiscard]] std::tuple<std::array<T,9>, std::array<R,3>>
    tuple() const {
        std::array<T,9> m{{xx, xy, xz, yx, yy, yz, zx, zy, zz}};
        std::array<R,3> v{{x, y, z}};
        return std::make_tuple(m, v);
    }
  };

  // Special simple version of Matrix for HighFive HDF5 interface
  template<class T>
  struct HF_Matrix {
      T xx, xy, xz, yx, yy, yz, zx, zy, zz;
  public:
      explicit HF_Matrix() = default;
      explicit HF_Matrix(const std::array<T, 9>& m):
      xx(m[0]), xy(m[1]), xz(m[2]), yx(m[3]), yy(m[4]), yz(m[5]), zx(m[6]), zy(m[7]), zz(m[8]) {}
      [[nodiscard]] std::array<T,9> array() const {
          return std::array<T,9>({xx, xy, xz, yx, yy, yz, zx, zy, zz});
      }
  };
} // namespace brille

// predeclare the specialisations of HighFive's create_datatype function
namespace HighFive {
template<> DataType create_datatype<brille::RotatesLike>();
template<> DataType create_datatype<brille::Bravais>();
template<> DataType create_datatype<brille::LengthUnit>();
template<> DataType create_datatype<brille::NodeType>();
template<> DataType create_datatype<brille::HF_Matrix<int>>();
template<> DataType create_datatype<brille::HF_Matrix<double>>();
template<> DataType create_datatype<brille::HF_Motion<int,double>>();
}

// templated functions for used with create_datatypes above, in case we want
// more than Matrix<int> and Motion<int,double> exposed to HighFive
template<class T, class R>
HighFive::CompoundType create_compound_Motion(){
  return {
      {"xx", ::HighFive::AtomicType<T>{}},
      {"xy", ::HighFive::AtomicType<T>{}},
      {"xz", ::HighFive::AtomicType<T>{}},
      {"yx", ::HighFive::AtomicType<T>{}},
      {"yy", ::HighFive::AtomicType<T>{}},
      {"yz", ::HighFive::AtomicType<T>{}},
      {"zx", ::HighFive::AtomicType<T>{}},
      {"zy", ::HighFive::AtomicType<T>{}},
      {"zz", ::HighFive::AtomicType<T>{}},
      {"x",  ::HighFive::AtomicType<R>{}},
      {"y",  ::HighFive::AtomicType<R>{}},
      {"z",  ::HighFive::AtomicType<R>{}},
  };
}
template<class T>
HighFive::CompoundType create_compound_Matrix(){
  return {
      {"xx", ::HighFive::AtomicType<T>{}},
      {"xy", ::HighFive::AtomicType<T>{}},
      {"xz", ::HighFive::AtomicType<T>{}},
      {"yx", ::HighFive::AtomicType<T>{}},
      {"yy", ::HighFive::AtomicType<T>{}},
      {"yz", ::HighFive::AtomicType<T>{}},
      {"zx", ::HighFive::AtomicType<T>{}},
      {"zy", ::HighFive::AtomicType<T>{}},
      {"zz", ::HighFive::AtomicType<T>{}},
  };
}
#endif //USE_HIGHFIVE
#endif //BRILLE_HDF_INTERFACE_HPP
