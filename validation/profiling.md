# Extra profiling output from `brille`

Profiling the Python/C++ interface in `brille` is challenging since neither a C++ or Python profiler can identify bottlenecks in the interface itself.
A Python profiler can not peer within the binary module and a C++ profiler can not examine the python code.
To help identify performance issues caused by, e.g., memory copies between the Python heap and C++ heap a special macro-bases profiling output is available in `brille`.
The profiling macros also print the system time in `YYYY-MM-DD HH:MM:SS:mmm` format by default when called, this can be useful for determining durations spanning Python and C++.

## Programmer's guide to `brille` output macros
A single logging object exists within the `brille` module which is used to provide output to the standard output through a set of macros representing different log-levels, `info`, `debug`, `verbose`, and `profile`.

For each log level there are two macros defined, `[level]_update(...)` and `[level]_update_if(logical_expression,...)` where `[level]` should be replaced by `info`, `debug`, `verbose`, or `profile`.
The first macro takes any number of arguments and uses them as the arguments for the logging object's print method.
The second macro checks whether its first argument evaluates to `true` and forwards the remaining arguments if so.

The macros are always defined but, other than the `info` macros, their definitions are empty unless a preprocessor variable matching their level is defined, i.e., `DEBUG`, `VERBOSE`, `PROFILE`; except that `VERBOSE` implies `DEBUG` as well.
When a macro is left as an empty definition its contents are elided away at compile-time, and the compiled executable does not execute any instructions associated with the macro.
A final macro `debug_exec(...)` exists purely to elide expensive operations in order to support more-complex debugging functionality which can be left in place once an operation is verified.

The log level can be controlled by the CMake configuration variable `BRILLE_LOGLEVEL`, e.g., `cmake .. -DBRILLE_LOGLEVEL=VERBOSE` will define the appropriate preprocessor variable at compile-time.
As profiling is potentially orthogonal to logging it is controlled by a separate CMake variable, `BRILLE_PROFILING`, and can be enabled via, e.g., `cmake .. -DBRILLE_PROFILING=On`

## Profiling output from the binary Python module
At present creating the `_brille` Python module with an elevated log level or profiling output requires *either* modifying `setup.py` to include the appropriate CMake variables (by adding to the `cmake_args` list) before running, e.g., `python setup.py develop`*or* configure CMake and build the module 'by hand', e.g., `cmake .. -DBRILLE_PROFILING=ON; cmake --build . --target _brille`, and then move the resulting module into the appropriate* location (*it's left as an exercise to determine what 'appropriate' means for your system).

# A profiling example
As originally written the `brille` module did not have a method by which to share memory between the C++ allocated heap and the Python allocated heap.
This led to unnecessary copying of data from Python into `brille` (eigenvalues, eigenvectors, Q points, etc.) and from `brille` back into Python (interpolated results).
This copying could represent a significant proportion of the time `brille` takes to interpolate at a provided point and needed to be addressed by sharing heap memory.
In order to ensure that the copy time could be quantified a simple script was constructed to use in conjunction with strategically placed profiling print statements.

```python
import psutil
import numpy as np
import brilleu as bu
from pathlib import Path

# Use a CASTEP binary file to build a brilleu object
castep_file = str(Path('C:\\', 'data', 'euphonic-benchmarking', 'lzo', 'castep', 'La2Zr2O7.castep_bin'))
breu = bu.BrillEu.from_castep(castep_file, max_volume=0.00001, parallel=True, sort=True, log={}, emit=True)

#(Q,E) box limits chosen to match limits from
# 	cut_sqw(lzo_hh2_qe.sqw, [-6,0.05,-5], [2,0.5,10]);
#
qemax = np.array([[-2.8, 3.6, -1,  10.2]])
qemin = np.array([[-3.5, 2.9, -1.3, 1.8]])

# The cut described above contains 192308 (Q,E) voxels
# the exact positions of which (probably) don't matter;
# for fair-comparisons we should seed the random number generator
rng = np.random.mtrand.RandomState(1515)
def genpoints(n):
	return rng.rand(n,4)*(qemax-qemin)+qemin
# We (probably) can't calculate all points in one go, so chunk the calculation
chunk_size = psutil.virtual_memory().free // breu.grid.bytes_per_point // 2

nq = 100000 # chosen to fit in one chunk on the development system
total = 0
opts = {'resfun':'gauss', 'param':[1,2.]}
sqw = np.ndarray(nq)
while total < nq:
	if total + chunk_size > nq:
		chunk_size = nq - total
	qe = genpoints(chunk_size)
	sqw[total:total+chunk_size] = breu.s_qw(qe[:,:-1], qe[:,-1], opts)
	total += chunk_size

print(breu.log)
```

When run with `brille == v0.4.3` this produced

```bash
[2020-10-30 17:45:41:778] Start of BrillouinZone construction
[2020-10-30 17:45:41:811]   End of BrillouinZone construction
[2020-10-30 17:45:41:811] Start of PolyhedronTrellis construction
[2020-10-30 17:45:43:385]   End of PolyhedronTrellis construction
[2020-10-30 17:45:43:385] Start of PolyhedronTrellis permutation key collection
[2020-10-30 17:45:43:422]   End of PolyhedronTrellis permutation key collection
[2020-10-30 17:45:52:190] Start of 'fill' operation with cost information
[2020-10-30 17:45:52:673]   End of 'fill' operation with cost information
[2020-10-30 17:45:52:673] Start of 'sort' operation
[2020-10-30 17:46:08:696]   End of 'sort' operation
[2020-10-30 17:46:08:728] Start of 'ir_interpolate_at' operation
[2020-10-30 17:46:08:732] Q array coppied from Python to C++ heap
[2020-10-30 17:46:08:732] Start BrillouinZoneTrellis3::ir_interpolate_at
[2020-10-30 17:46:08:733] BrillouinZone::ir_moveinto called with 6 threads
[2020-10-30 17:46:08:735] BrillouinZone::moveinto called with 6 threads
[2020-10-30 17:46:15:874] Parallel PolyhedronTrellis::interpolation_at 100000 points with 6 threads
[2020-10-30 17:46:25:003] Apply rotations/permutations to interpolated results
[2020-10-30 17:46:25:006] Start InnerInterpolationData::rip_real method
[2020-10-30 17:46:25:006] Start InnerInterpolationData::rip_gamma_complex method
[2020-10-30 17:46:28:296]   End InnerInterpolationData::rip_gamma_complex method
[2020-10-30 17:46:28:296]   End BrillouinZoneTrellis3::ir_interpolate_at
[2020-10-30 17:46:33:319] Interpolated values found
[2020-10-30 17:46:35:747]   End of 'ir_interpolate_at' operation
```

which captures three memory copies

1. the `fill` operation, 483 msec
2. copying Q points, 4 msec
3. returning the interpolated results, 7451 msec

Switching to shared memory arrays yielded

```bash
[2020-11-01 09:02:03:896] Start of BrillouinZone construction
[2020-11-01 09:02:03:920]   End of BrillouinZone construction
[2020-11-01 09:02:03:920] Start of PolyhedronTrellis construction
[2020-11-01 09:02:05:636]   End of PolyhedronTrellis Construction
[2020-11-01 09:02:05:636] Start of PolyhedronTrellis permutation key collection
[2020-11-01 09:02:05:669]   End of PolyhedronTrellis permutation key collection
[2020-11-01 09:02:14:294] Start of 'fill' operation with cost information
[2020-11-01 09:02:14:295]   End of 'fill' operation with cost information
[2020-11-01 09:02:14:295] Start of 'sort' operation
[2020-11-01 09:02:30:696]   End of 'sort' operation
[2020-11-01 09:02:30:713] Start of 'ir_interpolate_at' operation
[2020-11-01 09:02:30:713] Q array wrapped for C++ use
[2020-11-01 09:02:30:713] Start BrillouinZoneTrellis3::ir_interpolate_at
[2020-11-01 09:02:30:714] BrillouinZone::ir_moveinto called with 6 threads
[2020-11-01 09:02:30:720] BrillouinZone::moveinto called with 6 threads
[2020-11-01 09:02:49:413] Parallel interpolation at 100000 points with 6 threads
[2020-11-01 09:02:53:668] Apply rotations/permutations to interpolated results
[2020-11-01 09:02:53:673] Start Interpolator::rip_real method
[2020-11-01 09:02:53:673] Start Interpolator::rip_gamma_complex method
[2020-11-01 09:02:56:681]   End Interpolator::rip_gamma_complex method
[2020-11-01 09:02:56:681]   End BrillouinZoneTrellis3::ir_interpolate_at
[2020-11-01 09:02:56:682]   End of 'ir_interpolate_at' operation
```

with equivalent operations taking 1, 0, and 1 msec, respectively.
