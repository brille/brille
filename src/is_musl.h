/* brille needs a way to identify if it is being compiled against musl libc
 * since something in that library interferes with OpenMP and TetGen to cause
 * a segmentation fault in one brille test.
 *
 * Since the musllinux wheels are mainly used for documentation at the moment,
 * forcing OpenMP to use a single thread in case of musl libc in the offending
 * method will hopefully be sufficient for the time being
 * */
#ifndef BRILLE_IS_MUSL_H
#define BRILLE_IS_MUSL_H
/* from https://stackoverflow.com/a/70211227
 *
 * musl libc doesn't define __USE_GNU in features.h when _GNU_SOURCE is defined
 * and instead uses _GNU_SOURCE to guard GNU extensions.
 * So this procedure should work to create a musl libc macro
 * (it will have false positives if another libc implementation has the same
 * behavior or someone without a libc includes a stub features.h for some
 * reason; but Glibc, uClibc, and bionic all define __USE_GNU in features.h or
 * another header included by features.h when _GNU_SOURCE is defined)
 */
//#define _GNU_SOURCE 1
//#include <features.h>
//#ifndef __USE_GNU
//#define __MUSL__
//#endif
//#undef _GNU_SOURCE

#endif // BRILLE_IS_MUSL_H
