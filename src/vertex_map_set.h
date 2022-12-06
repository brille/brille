//
// Created by g on 30/11/22.
//

#ifndef BRILLE_VERTEX_MAP_SET_H
#define BRILLE_VERTEX_MAP_SET_H

#include "types.hpp"
#include "comparisons.hpp"
#include <vector>
#include <mutex>

namespace brille {

enum class AddVertexType {Unknown, Pristine, Crafted};

// Originally intended to be used with multi-threaded access
// Re-worked for each thread to have its own; so the mutexes would be
// superfluous and prevent move assignment of std::pair<VertexMapSet, R>

template <class T, template<class> class A>
class VertexMapSet {
public:
  using preserve_t = std::vector<ind_t>;
private:
  A<T> pristine_;
  A<T> appended_;
  preserve_t preserved_;
  ind_t appended_count_ = 0;
  ind_t preserved_count_ = 0;
//  std::mutex append_lock_;
//  std::mutex preserve_lock_;
  T relative_ = T(0);
  int digits_ = 0;

public:
  VertexMapSet(const A<T>& p, T s, int d): pristine_(p), relative_(s), digits_(d)
  {
    construct_preserved_appended();
  }
  VertexMapSet(A<T>&& p, T s, int d): pristine_(std::move(p)), relative_(s), digits_(d)
  {
    construct_preserved_appended();
  }
  VertexMapSet(A<T>&& p, preserve_t && pr, T r, int d):
    pristine_(std::move(p)), preserved_(std::move(pr)), relative_(r), digits_(d)
  {
    construct_appended();
  }
  VertexMapSet(A<T> && p, A<T> && a, preserve_t && pr, T r, int d):
    pristine_(std::move(p)), appended_(std::move(a)), preserved_(std::move(pr)),
    relative_(r), digits_(d)
  {
    appended_count_ = appended_.size(0);
    preserved_count_ = static_cast<ind_t>(
        std::count_if(preserved_.begin(), preserved_.end(),
                      [no=pristine_.size(0)](ind_t x){return x < no;}
                      ));
  }

  T relative() const {return relative_;}
  int digits() const {return digits_;}

  A<T> pristine() const {return pristine_;}
  A<T> origin() const {return T(0) * pristine_.view(0);}

  A<T> appended() const {return appended_.view(0, appended_count_);}

  [[nodiscard]] ind_t pristine_count() const { return pristine_.size(0);}
  [[nodiscard]] ind_t preserved_count() const{ return preserved_count_;}
  [[nodiscard]] ind_t appended_count() const { return appended_count_; }

  bool reserve_appended(ind_t total) {
    if (total > appended_.size(0)){
      appended_.resize(total);
      return true;
    }
    return false;
  }

  [[nodiscard]] std::tuple<bool, ind_t> pristine_origin() const {
    return is_pristine(origin());
  }
  [[nodiscard]] bool is_preserved(ind_t i) const {return preserved_[i] < pristine_.size(0);}
  [[nodiscard]] ind_t preserved_value(ind_t i) const {return preserved_[i];}
  [[nodiscard]] const preserve_t& preserved() const {return preserved_;}

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
    auto count = equals.count(appended_count_); // avoid matching uninitialized entries beyond appended_count_
    if (count > 1) throw std::runtime_error("Too many extra matches to vertex!");
    return count ? std::make_tuple(true, pristine_count() + equals.first(appended_count_)) : std::make_tuple(false, ind_t(0));
  }

  ind_t preserve(ind_t index) {
//    std::lock_guard<std::mutex> guard(preserve_lock_);
    if (preserved_[index] > pristine_count()) preserved_[index] = preserved_count_++;
    return preserved_[index];
  }

  ind_t add(const A<T>& vertex, AddVertexType type = AddVertexType::Unknown) {
    if (type != AddVertexType::Crafted){
      auto [pristine_present, pristine_index] = is_pristine(vertex);
      if (pristine_present) return preserve(pristine_index);
    }
    // If the user *thought* the point was pristine, but it was not
    // it will be added to here ... Maybe not great.
    auto [present, index] = is_crafted(vertex);
    return present ? index : insert(vertex);
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

  ind_t origin_index() {
    // Return the index of the origin, add it to appended_ if necessary
    return add(origin());
  }

  VertexMapSet<T,A> consolidate() const {
    std::vector<size_t> extract;
    extract.reserve(preserved_count_);
    for (size_t i = 0; i < preserved_count_; ++i){
      auto itr = std::find(preserved_.begin(), preserved_.end(), i);
      debug_update_if(itr == preserved_.end(), "Could not find index with preserved value ", i, "?!");
      if (itr != preserved_.end())
      extract.push_back(std::distance(preserved_.begin(), itr));
      else
      throw std::runtime_error("Could not find the value!");
    }
    // preserved_ [0, 3, 1, x, x, 2, 4, 5, …] extracts [0, 2, 5, 1, 6, 7, …]
    // from pristine_;
    // combine the retained vertices and the appended vertices,
    // protecting against view(0, 0)
    if (appended_count_) {
      auto nv = cat(0, pristine_.extract(extract), appended());
      preserve_t pr;
      auto no = nv.size(0);
      pr.resize(no);
      std::iota(pr.begin(), pr.end(), 0);
      return {std::move(nv), std::move(pr), relative_, digits_};
    } else {
      auto nv = pristine_.extract(extract);
      preserve_t pr;
      auto no = nv.size(0);
      pr.resize(no);
      std::iota(pr.begin(), pr.end(), 0);
      return {std::move(nv), std::move(pr), relative_, digits_};
    }
  }

  [[nodiscard]] bool is_consolidated() const {
    if (appended_.size(0) > 0) return false;
    for (ind_t i=0; i<preserved_.size(); ++i) if (preserved_[i] != i) return false;
    return true;
  }

  friend std::ostream & operator<<(std::ostream & os, const VertexMapSet<T,A>& v){
    os << "VertexMapSet\n";
    os << "\tpristine: (" << v.pristine().size(0) << "->" << v.preserved_count() << ")\n";
    os << "\tappended: (" << v.appended_count() << ")\n";
    os << "\tpreserved:\n";
    for (const auto & p: v.preserved()) os << p << " ";
    os << "\n";
    return os;
  }

private:
  void construct_preserved_appended() {
    preserved_ = std::vector<ind_t>(pristine_.size(0), pristine_.size(0)+1u);
    construct_appended();
  }
  void construct_appended() {
    appended_ = 0 * pristine_.view(0); // 0 * forces a new array to be created
    appended_.resize(pristine_.size(0) >> 1u);
  }

  ind_t insert(const A<T>& vertex){
//    std::lock_guard<std::mutex> guard(append_lock_);
    if (appended_.size(0) < appended_count_ + 1u) appended_.resize(2 * appended_count_);
    appended_.set(appended_count_, vertex);
    return pristine_count() + appended_count_++;
  }
};

} // namespace brille

#endif // BRILLE_VERTEX_MAP_SET_H
