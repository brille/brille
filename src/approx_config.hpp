#ifndef BRILLE_APPROX_CONFIG_HPP_
#define BRILLE_APPROX_CONFIG_HPP_

#include <cmath>
#include <memory>
namespace brille::approx_float {
  class ApproxConfig {
    int digit_; //! default 1000
    int power_; //! default -15; the floating point tolerance in powers of 10
  public:
    explicit ApproxConfig(): digit_{1000}, power_{-15} {}
    ApproxConfig(int dig, int pow): digit_{dig}, power_{pow} {}
    [[nodiscard]] int digit() const {return digit_;};
    [[nodiscard]] int power() const {return power_;}
    template<class T> [[nodiscard]] T tol() const {return std::pow(T(10), static_cast<T>(power_));}
    int digit(int i) {digit_=i; return digit_;}
    int power(int i) {power_=i; return power_;}
  };
}
//namespace brille {
//[[maybe_unused]] auto approx_config = std::make_shared<approx_float::ApproxConfig>();
//}

#endif