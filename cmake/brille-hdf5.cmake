# HighFive is used to simplify the interface with HDF5:

# Conan handles getting HighFive into the right place with the right build options:
#   No support for Boost, Eigen, or Xtensor
find_package(HighFive REQUIRED)
foreach(HF_TARGET IN LISTS CXX_TARGETS)
    target_link_libraries(${HF_TARGET} PUBLIC HighFive)
endforeach()