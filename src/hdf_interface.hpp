//
// Created by g on 19/11/2021.
//

#ifndef BRILLE_HDF_INTERFACE_HPP
#define BRILLE_HDF_INTERFACE_HPP

#include <map>
#include <vector>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

const std::string LENGTH_ENTRY("/length");
const std::string LENGTH_ATTR_NAME("length");

namespace brille {
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
  lists_from_hdf(const R& obj, const std::string& entry){
    using namespace HighFive;
    auto group = obj.getGroup(entry);
    auto attr = group.getAttribute(LENGTH_ATTR_NAME);
    size_t no;
    attr.read(no);
    std::vector<std::vector<T>> out(no);
    for (size_t i=0; i<no; ++i){
        DataSet dataset = group.getDataSet(std::to_string(i));
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
    std::vector<Key> keys(map.size());
    std::vector<Value> values(map.size());
    for (const auto [key, value]: map){
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

}
#endif //BRILLE_HDF_INTERFACE_HPP
