# The hdf5 recipe builds with cmake/[>=3.18 <4], which has no generator for
# Visual Studio 2026 (msvc 195); that needs CMake 4.2. A recipe cannot override
# a dependency's tool_requires, but a profile can.
[replace_tool_requires]
cmake/*: cmake/[>=4.2 <5]
