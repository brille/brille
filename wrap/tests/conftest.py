"""Pin BLAS to one thread before numpy is imported.

The reference models diagonalise many small matrices, where a threaded BLAS
only adds overhead. Worse, OpenBLAS's pthreads pool spin-waits and competes
with brille's OpenMP threads: on a 6-core machine one test took over 300 s
instead of 5 s. Set the variables yourself to override.
"""
import os

for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ.setdefault(_name, "1")
