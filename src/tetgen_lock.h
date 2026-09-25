#ifndef BRILLE_TETGEN_LOCK_H_
#define BRILLE_TETGEN_LOCK_H_
#include <mutex>

namespace brille {
  /*! \brief The lock every call into TetGen must hold

  The bundled TetGen is not thread safe: tetrahedralizing from several threads at
  once gives different (or broken) meshes. With the GIL released, Python threads
  can build meshes concurrently, so all calls into TetGen are serialized here.
  */
  inline std::mutex & tetgen_mutex() {
    static std::mutex mutex;
    return mutex;
  }
}
#endif
