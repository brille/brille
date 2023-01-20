//
// Created by g on 30/11/22.
//

#ifndef BRILLE_VERTEX_MAP_SET_H
#define BRILLE_VERTEX_MAP_SET_H

#include "types.hpp"
#include "comparisons.hpp"
#include <vector>
#include <mutex>
#include <map>

namespace brille {

enum class AddVertexType {Unknown, Pristine, Crafted};
enum class MapVertexType {Pristine, Appended, SecondPristine, SecondAppended};

// Originally intended to be used with multi-threaded access
// Re-worked for each thread to have its own; so the mutexes would be
// superfluous and prevent move assignment of std::pair<VertexMapSet, R>

template <class T, template<class> class A>
class VertexMapSet {
public:
//  using preserve_t = std::vector<ind_t>;
  using rel_t = std::pair<MapVertexType, ind_t>;
  using preserve_t = std::map<ind_t, ind_t>;
private:
  A<T> pristine_;
  A<T> appended_;
  preserve_t preserve_;
  ind_t appended_count_ = 0;
  ind_t preserve_count_ = 0;
//  std::mutex append_lock_;
//  std::mutex preserve_lock_;
  T relative_ = T(0);
  int digits_ = 0;

public:
  explicit VertexMapSet() = default;
  VertexMapSet(const A<T>& p, T s, int d): pristine_(p), relative_(s), digits_(d)
  {
    construct_preserved_appended();
  }
  VertexMapSet(A<T>&& p, T s, int d): pristine_(std::move(p)), relative_(s), digits_(d)
  {
    construct_preserved_appended();
  }
  VertexMapSet(A<T>&& p, preserve_t && pr, T r, int d):
    pristine_(std::move(p)), preserve_(std::move(pr)), relative_(r), digits_(d)
  {
    construct_appended();
  }
  VertexMapSet(A<T> && p, A<T> && a, preserve_t && pr, T r, int d):
    pristine_(std::move(p)), appended_(std::move(a)), preserve_(std::move(pr)),
    relative_(r), digits_(d)
  {
    appended_count_ = appended_.size(0);
    preserve_count_ = static_cast<ind_t>(preserve_.size());
//    preserved_count_ = static_cast<ind_t>(
//        std::count_if(preserved_.begin(), preserved_.end(),
//                      [no=pristine_.size(0)](ind_t x){return x < no;}
//                      ));
  }

  T relative() const {return relative_;}
  [[nodiscard]] int digits() const {return digits_;}

  A<T> pristine() const {return pristine_;}
  A<T> origin() const {return T(0) * pristine_.view(0);}

  A<T> appended() const {return appended_.view(0, appended_count_);}

  [[nodiscard]] ind_t pristine_count() const { return pristine_.size(0);}
  [[nodiscard]] ind_t preserved_count() const{ return static_cast<ind_t>(preserve_.size());}
  [[nodiscard]] ind_t appended_count() const { return appended_count_; }

  [[nodiscard]] const preserve_t& preserved() const {return preserve_;}

  void reserve_appended(ind_t total) {
    if (total > appended_.size(0)) appended_.resize(total);
  }

private:
  [[nodiscard]] std::tuple<bool, ind_t> pristine_origin() const {
    return is_pristine(origin());
  }
  [[nodiscard]] bool is_preserved(ind_t i) const {
    const auto & p{preserve_}; // force the use of const_iterator
    auto search = p.find(i);
    return search != p.end();
  }
  [[nodiscard]] ind_t preserved_value(ind_t i) const {
    const auto & p{preserve_}; // force the use of const_iterator
    if(auto search = p.find(i); search != p.end())
      return search->second;
    throw std::runtime_error("Out-of-bounds preserved valued");
  }
//  [[nodiscard]] bool is_preserved(ind_t i) const {return preserved_[i] < pristine_.size(0);}
//  [[nodiscard]] ind_t preserved_value(ind_t i) const {return preserved_[i];}

  std::tuple<bool, ind_t> is_pristine(const A<T>& vertex) const {
    // pristine can not change; so no need to lock its mutex:
    auto equals = pristine_.row_is(brille::cmp::eq, vertex, relative_, relative_, digits_);
    auto count = equals.count();
    if (count > 1) throw std::runtime_error("Too many matches to vertex!");
    return count ? std::make_tuple(true, equals.first()) : std::make_tuple(false, ind_t(0));
  }

  std::tuple<bool, ind_t> is_crafted(const A<T>& vertex) {
    // appended and appended_count *CAN CHANGE* so we must lock the mutex
//    std::lock_guard<std::mutex> guard(append_lock_);

    auto equals = appended_.row_is(brille::cmp::eq, vertex, relative_, relative_, digits_);
    auto count = (appended_count_ > 0 ? equals.count(appended_count_) : 0); // avoid matching uninitialized entries beyond appended_count_
    if (count > 1) {
      std::cout << vertex << " has " << count << " matches in first " << appended_count_ << " entries of \n" << appended_ << " but only 1 expected" << std::endl;
      throw std::runtime_error("Too many extra matches to vertex!");
    }
//    return count ? std::make_tuple(true, pristine_count() + equals.first(appended_count_)) : std::make_tuple(false, ind_t(0));
    return count ? std::make_tuple(true, equals.first(appended_count_)) : std::make_tuple(false, ind_t(0));
  }

public:
  rel_t preserve(ind_t index) {
    if (!is_preserved(index)) preserve_[index] = static_cast<ind_t>(preserve_.size());
    return std::make_pair(MapVertexType::Pristine, preserve_[index]);
  }

  rel_t add(const A<T>& vertex, AddVertexType add_type = AddVertexType::Unknown) {
    if (add_type != AddVertexType::Crafted){
      auto [pristine_present, pristine_index] = is_pristine(vertex);
      if (pristine_present) return preserve(pristine_index);
    }
    // If the user *thought* the point was pristine, but it was not
    // it will be added to here ... Maybe not great.
    auto [present, index] = is_crafted(vertex);
    return present ? std::make_pair(MapVertexType::Appended, index) : insert(vertex);
  }

//  ind_t add(A<T>&& vertex, AddVertexType type = AddVertexType::Unknown){
//    if (type != AddVertexType::Crafted){
//      auto [pristine_present, pristine_index] = is_pristine(vertex);
//      if (pristine_present) return preserve(pristine_index);
//    }
//    // If the user *thought* the point was pristine, but it was not
//    // it will be added to here ... Maybe not great.
//    auto [present, index] = is_crafted(vertex);
//    return present ? index : insert(vertex);
//  }

  rel_t origin_index() {
    // Return the index of the origin, add it to appended_ if necessary
    return add(origin());
  }



  VertexMapSet<T,A> consolidate() const {
    std::vector<size_t> extract;
    extract.reserve(preserve_.size());
    // Pull out vertices based on the order in which they were preserved
    for (const auto & pair: key_sorted_preserve()){
      extract.push_back(pair.first);
    }
//    for (size_t i = 0; i < preserved_count_; ++i){
//      auto itr = std::find(preserved_.begin(), preserved_.end(), i);
//      debug_update_if(itr == preserved_.end(), "Could not find index with preserved value ", i, "?!");
//      if (itr != preserved_.end())
//      extract.push_back(std::distance(preserved_.begin(), itr));
//      else
//      throw std::runtime_error("Could not find the value!");
//    }
    // preserve_ {{0,0}, {1,3}, {2,1}, {5,2}, {6,4}, {7,5}, …}
    // is sorted to [{0,0}, {2,1}, {5,2}, {1,3}, {6,4}, {7,5}, …]
    // which then extracts [0, 2, 5, 1, 6, 7, …] from pristine_

    // preserved_ [0, 3, 1, x, x, 2, 4, 5, …] extracts [0, 2, 5, 1, 6, 7, …]
    // from pristine_;
    // combine the retained vertices and the appended vertices,
    // protecting against view(0, 0)
    if (appended_count_) {
      auto nv = cat(0, pristine_.extract(extract), appended());
      preserve_t pr;
      auto no = nv.size(0);
      for (ind_t i=0; i<no; ++i) pr[i] = i;
      return {std::move(nv), std::move(pr), relative_, digits_};
    } else {
      auto nv = pristine_.extract(extract);
      preserve_t pr;
      auto no = nv.size(0);
      for (ind_t i=0; i<no; ++i) pr[i] = i;
      return {std::move(nv), std::move(pr), relative_, digits_};
    }
  }

  [[nodiscard]] bool is_consolidated() const {
    if (appended_.size(0) > 0) return false;
//    for (ind_t i=0; i<preserved_.size(); ++i) if (preserved_[i] != i) return false;
    for (const auto & pr: preserve_) if (pr.first != pr.second) return false;
    return true;
  }

  friend std::ostream & operator<<(std::ostream & os, const VertexMapSet<T,A>& v){
    os << "VertexMapSet\n";
    os << "  pristine: (" << v.pristine().size(0) << "->" << v.preserved_count() << ")\n";
//    os << v.pristine();
    os << "  appended: (" << v.appended_count() << ")\n";
//    os << v.appended();
//    os << "\tpreserved:\n    ";
    os << "  preserved: ";
    for (const auto & p: v.preserved()) os << "{" << p.first << ":" << p.second << "} ";
//    os << "\n";
    return os;
  }

private:
  void construct_preserved_appended() {
//    preserved_ = std::vector<ind_t>(pristine_.size(0), pristine_.size(0)+1u);
    construct_appended();
  }
  void construct_appended() {
    appended_ = 0 * pristine_.view(0); // 0 * forces a new array to be created
    appended_.resize(pristine_.size(0) >> 1u);
  }

  rel_t insert(const A<T>& vertex){
    if (appended_.size(0) < appended_count_ + 1u) appended_.resize(2 * appended_count_);
    appended_.set(appended_count_, vertex);
    //return pristine_count() + appended_count_++;
    return std::make_pair(MapVertexType::Appended, appended_count_++);
  }

  [[nodiscard]] std::vector<std::pair<ind_t, ind_t>> key_sorted_preserve() const {
    std::vector<std::pair<ind_t, ind_t>> ksp;
    for (auto & pair: preserve_) ksp.push_back(pair);
    std::sort(ksp.begin(), ksp.end(), [](auto & a, auto & b){return a.second < b.second;});
    return ksp;
  }
};
} // namespace brille

std::ostream & operator<<(std::ostream & os, brille::MapVertexType t);

#endif // BRILLE_VERTEX_MAP_SET_H
