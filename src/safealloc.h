/*! \file */
#ifndef _SAFEALLOC_H_
#define _SAFEALLOC_H_
// (array) allocation for specified type

/*! \brief safer(?) array allocation for a specified type */
template<typename R> R * safealloc(const size_t d){
  R *out = nullptr;
  out = new R[d](); // R[d] does not initialize to zero. R[d]() does!
  if (out == nullptr) throw std::bad_alloc();
  return out;
}
#endif
