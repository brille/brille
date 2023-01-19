set(CONAN_LLVM_OPENMP "")
if (APPLE)
    set(CONAN_LLVM_OPENMP llvm-openmp/11.1.0)
endif()

# Use Conan to fetch/build HDF5, the manylinux2014 image has HDF5 v1.8 which causes errors in testing
if (NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
            "${CMAKE_BINARY_DIR}/conan.cmake"
            EXPECTED_HASH SHA256=5cdb3042632da3efff558924eecefd580a0e786863a857ca097c3d1d43df5dcd
            TLS_VERIFY ON)
endif()
include("${CMAKE_BINARY_DIR}/conan.cmake")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_BINARY_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})
if (DEFINED ENV{CONAN_USER_HOME})
    message(STATUS "Use $ENV{CONAN_USER_HOME} as CONAN_USER_HOME")
endif()
conan_cmake_configure(
        REQUIRES
        hdf5/1.12.0
        ${CONAN_LLVM_OPENMP}
        GENERATORS
        cmake_find_package
        OPTIONS
        hdf5:shared=False
        hdf5:hl=False
        hdf5:with_zlib=False
)
conan_cmake_autodetect(conan_settings)
conan_cmake_install(PATH_OR_REFERENCE ${CMAKE_CURRENT_BINARY_DIR} BUILD outdated REMOTE conancenter SETTINGS ${conan_settings})
