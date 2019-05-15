#ifndef _UNSIGNEDTOSIGNED_H
#define _UNSIGNEDTOSIGNED_H

#include <limits>     // for std::numeric_limits
#include <stdexcept>  // for std::overflow_error
#include <string>

template<typename S,typename U>
S unsigned_to_signed(const U u){
  if (u > static_cast<U>((std::numeric_limits<S>::max)())){
    std::string msg = "Value " + std::to_string(u)
                    + " can not be stored in requested signed type"
                    + " (max value "+ std::to_string((std::numeric_limits<S>::max)()) + ")";
    throw std::overflow_error(msg);
  }
  return static_cast<S>(u);
}

template<typename U,typename S>
U signed_to_unsigned(const S s){
  if (s < 0 || static_cast<S>((std::numeric_limits<U>::max)())){
    std::string msg = "Value " + std::to_string(s)
                    + " can not be stored in requested unsiged type"
                    + " (max value "+ std::to_string((std::numeric_limits<U>::max)()) + ")";
  }
  return static_cast<U>(s);
}


#endif
