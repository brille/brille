#ifndef BRILLE_BZ_CONFIG_HPP_
#define BRILLE_BZ_CONFIG_HPP_

#include <cmath>
namespace brille {
  class BrillouinZoneConfig {
    bool primitive_; //! default true
    bool time_reversal_; //! default false
    bool wedge_search_; //! default true
    bool divide_primitive_; //! default true
    int divide_extent_; //! default 1 (will use 2 as well)
  public:
    explicit BrillouinZoneConfig():
      primitive_{true}, time_reversal_{false}, wedge_search_{true},
      divide_primitive_{true}, divide_extent_{1} {}
    [[nodiscard]] bool primitive() const {return primitive_;}
    [[nodiscard]] bool time_reversal() const {return time_reversal_;}
    [[nodiscard]] bool wedge_search() const {return wedge_search_;}
    [[nodiscard]] bool divide_primitive() const {return divide_primitive_;}
    [[nodiscard]] int divide_extent() const {return divide_extent_;}
    bool primitive(bool p) {primitive_ = p; return primitive_;}
    bool time_reversal(bool t) {time_reversal_ = t; return time_reversal_;}
    bool wedge_search(bool p) {wedge_search_ = p; return wedge_search_;}
    bool divide_primitive(bool p) {divide_primitive_ = p; return divide_primitive_;}
    int divide_extent(int i) {divide_extent_=i; return divide_extent_;}
  };
}

#endif