#ifndef BRILLE_APPROX_CONFIG_HPP_
#define BRILLE_APPROX_CONFIG_HPP_

#include <cstdint>
#include <cmath>
#include <memory>
#include "hdf_interface.hpp"

namespace brille::approx_float {
class Config {
private:
  int digit_;
  double direct_;
  double reciprocal_;

public:
  Config(): digit_{1000}, direct_{1e-12}, reciprocal_{1e-12} {}
  Config(int dig, double dir, double rec): digit_{dig}, direct_{dir}, reciprocal_{rec} {}

  // getters
  [[nodiscard]] int digit() const {return digit_;}
  template<class T> [[nodiscard]] T direct() const {return static_cast<T>(direct_);}
  template<class T> [[nodiscard]] T reciprocal() const {return static_cast<T>(reciprocal_);}

  // setters returning a set object reference (to, hopefully, allow chaining)
  Config & digit(int i) {digit_=i; return *this;}
  template<class T> Config direct(T value){
    direct_ = static_cast<double>(value);
    return *this;
  }
  template<class T> Config & reciprocal(T value){
    reciprocal_ = static_cast<double>(value);
    return *this;
  }
  template<class T> Config & tol(T value) {
    direct_ = static_cast<double>(value);
    reciprocal_ = direct_;
    return *this;
  }


#ifdef USE_HIGHFIVE
  template<class H> std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool>
  to_hdf(H& obj, const std::string & entry) const {
    auto group = overwrite_group(obj, entry);
    group.createAttribute("digit", digit_);
    group.createAttribute("direct", direct_);
    group.createAttribute("reciprocal", reciprocal_);
    return true;
  }
  template<class H> static std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, Config>
  from_hdf(H& obj, const std::string & entry) {
    int d{0};
    double dir{0}, rec{0};
    auto group = obj.getGroup(entry);
    group.getAttribute("digit").read(d);
    group.getAttribute("direct").read(dir);
    group.getAttribute("reciprocal").read(rec);
    return {d, dir, rec};
  }
#endif // USE_HIGHFIVE
};
//
//class ConfigHolder {
//public:
//  static Config& get();
//private:
//  ConfigHolder() = default;
//  ~ConfigHolder() = default;
//  ConfigHolder(ConfigHolder&) = delete;
//  ConfigHolder& operator=(ConfigHolder&) = delete;
//};

extern Config config;
}

#endif