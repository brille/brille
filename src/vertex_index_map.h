//
// Created by g on 02/12/22.
//

#ifndef BRILLE_VERTEX_INDEX_MAP_H
#define BRILLE_VERTEX_INDEX_MAP_H

#include "vertex_map_set.h"

namespace brille {

class VertexIndexMap{
public:
  using rel_t = std::pair<MapVertexType, ind_t>;
  using map_t = std::vector<rel_t>;
  using data_t = std::map<ind_t, map_t>;
protected:
  data_t data_;

public:
  explicit VertexIndexMap(): data_() {};
  explicit VertexIndexMap(const data_t & d): data_() {
    for (const auto & pair: d){
      map_t v;
      v.reserve(pair.second.size());
      for (const auto & x: pair.second) v.emplace_back(x);
      data_[pair.first] = v;
    }
  }
  VertexIndexMap(data_t && d): data_(std::move(d)) {}

  VertexIndexMap first_copy() const {
    return VertexIndexMap(data_);
  }

  VertexIndexMap second_copy() const {
    data_t sd;
    for (const auto & pair: data_){
      map_t v;
      v.reserve(pair.second.size());
      for (const auto & x: pair.second){
        MapVertexType mvt;
        switch(x.first){
        case MapVertexType::Pristine: mvt = MapVertexType::SecondPristine; break;
        case MapVertexType::Appended: mvt = MapVertexType::SecondAppended; break;
        default: throw std::runtime_error("Second copy of a second copy not allowed");
        }
        v.emplace_back(mvt, x.second);
      }
      sd[pair.first] = v;
    }
    return {std::move(sd)};
  }

  [[nodiscard]] bool has_second() const {
    for (const auto & pair: data_) for (const auto & x: pair.second){
      if (x.first == MapVertexType::SecondAppended || x.first == MapVertexType::SecondPristine) return true;
    }
    return false;
  }

  [[nodiscard]] bool merge(VertexIndexMap& other){
    data_.merge(other.data_);
    return other.size() > 0;
  }

  [[nodiscard]] const data_t & data() const {return data_;}

  void append(ind_t node, rel_t value) { data_[node].push_back(value); }
  void set(ind_t node, ind_t vertex, rel_t value) { data_[node][vertex] = value; }

  [[nodiscard]] size_t size() const {return data_.size();}

//  map_t & init(ind_t node, ind_t len){
//    // Why does this cause a fatal error when multiple threads are used?
//      auto & v = data_[node]; // Constructs a map_t if it doesn't exist
//      v.resize(len);
//      return v;
//  }
  map_t & get(ind_t node) {return data_[node]; }

  [[nodiscard]] const map_t & get(ind_t node) const {return data_.at(node);}
  [[nodiscard]] rel_t get(ind_t node, ind_t vertex) const {
    const auto & d{data_};
    if (auto search=d.find(node); search != d.end()){
      if (vertex < search->second.size())
        return search->second[vertex];
      throw std::runtime_error("Out-of-bounds vertex index!");
    }
    throw std::runtime_error("Node index not found!");
  }
  [[nodiscard]] ind_t decode(ind_t appended_offset, ind_t node, ind_t vertex) const {
    auto type_index = get(node, vertex);
    switch(type_index.first){
    case MapVertexType::Pristine: return type_index.second;
    case MapVertexType::Appended: return appended_offset + type_index.second;
    default: throw std::runtime_error("Only pristine and appended can be decoded");
    }
  }
  [[nodiscard]] std::vector<ind_t> decode(ind_t appended_offset, ind_t node) const {
    std::vector<ind_t> out;
    out.reserve(data_.at(node).size());
    for (const auto & x: data_.at(node)){
      switch(x.first){
      case MapVertexType::Pristine: out.push_back(x.second); break;
      case MapVertexType::Appended: out.push_back(appended_offset + x.second); break;
      default:  throw std::runtime_error("Only pristine and appended can be decoded");
      }
    }
    return out;
  }
  [[nodiscard]] bool has(ind_t node) const {
    const auto & d{data_};
    auto search=d.find(node);
    return search != d.end();
  }

  void replace(rel_t was, rel_t to){
    for (auto & pair: data_) {
      for (auto & val: pair.second) {
        if (val == was) {
          val = to;
        }
      }
    }
  }

  [[nodiscard]] typename data_t::const_iterator cbegin() const {return data_.cbegin();}
  [[nodiscard]] typename data_t::const_iterator cend() const {return data_.cend();}

  friend std::ostream & operator<<(std::ostream & os, const VertexIndexMap& v){
    for (auto it = v.cbegin(); it != v.cend(); ++it){
      os << it->first << ": [ ";
      for (auto vit = it->second.cbegin(); vit != it->second.cend(); ++vit){
        //os << vit->first << vit->second << " ";
				os << vit->second << " ";
			}
      os << "]\n";
    }
    os << "\n";
    return os;
  }
};

}

#endif // BRILLE_VERTEX_INDEX_MAP_H
