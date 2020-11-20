/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */
#ifndef BRILLE_SORTING_STATUS_HPP_
#define BRILLE_SORTING_STATUS_HPP_
namespace brille {
/*! \brief A type to hold sorting status information for an object

Every vertex/node/etc. in a grid-like object will have sorting information
associated with it. Notably, whether or not it is considered sorted, if it is
allowed to be further modified, and how many times its sorting permutation has
been changed.
The first two are flags and the last a positive integer, so this information
could all be stored in a single unsigned N-bit integer with the visit count
occupying only N-2 bits. This is probably overkill
*/
class SortingStatus{
public:
  using status_t = unsigned int;
private:
  status_t status;
  static const status_t s_mask{1u};
  static const status_t l_mask{2u};
  static const status_t q_mask{4u};
  static const status_t c_mask{8u};
  static const status_t v_mask{static_cast<status_t>(-1)^15u};
  static const status_t v_shift{4u};
public:
  explicit SortingStatus(status_t s=0u): status{s} {};
  SortingStatus(bool s, bool f, unsigned int v): status{0} {
    this->sorted(s);
    this->locked(f);
    this->visits(v);
  }
  bool sorted() const { return status & s_mask; }
  bool locked() const { return status & l_mask; }
  bool queued() const { return status & q_mask; }
  bool claimed() const { return status & c_mask; }
  bool visited() const { return (status & v_mask) > 0; }
  status_t visits() const { return status>>v_shift; }
  bool sorted(bool s) {
    if (this->sorted() != s) {
      status_t rem = this->status & (l_mask + q_mask + c_mask + v_mask);
      this->status = rem + (s ? s_mask : 0u);
    }
    return this->sorted();
  }
  bool locked(bool l){
    if (this->locked() != l){
      status_t rem = this->status & (s_mask + q_mask + c_mask + v_mask);
      this->status = rem + (l ? l_mask : 0u);
    }
    return this->locked();
  }
  bool queued(bool q){
    if (this->queued() != q){
      status_t rem = this->status & (s_mask + l_mask + c_mask + v_mask);
      this->status = rem + (q ? q_mask : 0u);
    }
    return this->queued();
  }
  bool claimed(bool c){
    if (this->claimed() != c){
      status_t rem = this->status & (s_mask + l_mask + q_mask + v_mask);
      this->status = rem + (c ? c_mask : 0u);
    }
    return this->claimed();
  }
  template<typename T>
  status_t visits(T v) {
    if (v <= static_cast<status_t>(-1) >> v_shift){
      this->status ^= v_mask; // keep everything but the number of visits
      this->status += v << v_shift; // add the offset visits value
    }
    else
      throw std::overflow_error("SortingStatus can not hold more than 2^29-1 visits.");
    return this->visits();
  }
  status_t addvisit() {
    if (this->status < v_mask)
      this->status += (1u << v_shift);
    else
      throw std::overflow_error("SortingStatus can not hold more than 2^29-1 visits.");
    return this->visits();
  }
  template<class T>
  status_t addvisit(const T max_visits){
    if (this->visits() < max_visits) return this->addvisit();
    this->locked(true);
    return this->visits();
  }
  template<class T>
  bool unlocked_addvisit_unsorted(const T max_visits){
    if (this->locked()) return false;
    this->addvisit(max_visits);
    if (this->sorted()) return false;
    return true;
  }
  std::string to_string() const {
    std::string str;
    bool s{this->sorted()}, l{this->locked()};
    if (s) str += "sorted";
    if (l){
      if (s) str += " and ";
      str += "locked";
    }
    if (l||s) str += " in ";
    str += std::to_string(this->visits()) + " visits";
    return str;
  }

  bool is_queable() const {
    return !(this->locked() || this->sorted() || this->queued() || this->claimed());
    // only queable if not locked, sorted, queued, or claimed
    // 15u == 1u + 2u + 4u + 8u
    // return 15u == (this->status&15u)^15u;
  }

};

}
#endif
