//
// Created by g on 02/12/22.
//

#ifndef BRILLE_VERTEX_INDEX_MAP_H
#define BRILLE_VERTEX_INDEX_MAP_H

namespace brille {

class VertexIndexMap{
public:
  using map_t = std::vector<ind_t>;
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

  bool merge(VertexIndexMap& other){
    data_.merge(other.data_);
    return other.size() > 0;
  }

  [[nodiscard]] const data_t & data() const {return data_;}

  void append(ind_t node, ind_t value) { data_[node].push_back(value); }
  void set(ind_t node, ind_t vertex, ind_t value) { data_[node][vertex] = value; }

  [[nodiscard]] size_t size() const {return data_.size();}

//  map_t & init(ind_t node, ind_t len){
//    // Why does this cause a fatal error when multiple threads are used?
//      auto & v = data_[node]; // Constructs a map_t if it doesn't exist
//      v.resize(len);
//      return v;
//  }
  map_t & get(ind_t node) {return data_[node]; }

  [[nodiscard]] const map_t & get(ind_t node) const {return data_.at(node);}
  [[nodiscard]] ind_t get(ind_t node, ind_t vertex) const {
    const auto & d{data_};
    if (auto search=d.find(node); search != d.end()){
      if (vertex < search->second.size())
        return search->second[vertex];
      throw std::runtime_error("Out-of-bounds vertex index!");
    }
    throw std::runtime_error("Node index not found!");
  }
  [[nodiscard]] bool has(ind_t node) const {
    const auto & d{data_};
    auto search=d.find(node);
    return search != d.end();
  }

  void replace(ind_t was, ind_t to){
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
      for (auto vit = it->second.cbegin(); vit != it->second.cend(); ++vit)
        os << *vit << " ";
      os << "]\n";
    }
    os << "\n";
    return os;
  }
};

}

#endif // BRILLE_VERTEX_INDEX_MAP_H
