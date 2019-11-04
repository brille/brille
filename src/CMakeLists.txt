# Directory specific compilation flags
if (MSVC)
    # warning level 4 -- add /WX for all warnings as errors
    # /Zc__cplusplus ensures the preprocessor macro __cplusplus is set
    add_compile_options(/W4 /Zc__cplusplus)
else()
    # lots of warnings -- add -Werror for  all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CXX_TARGETS _symbz test_cpp test_catch)

set(CXX_SOURCES
  arithmetic.cpp
  #balltree.cpp
  balltrellis.cpp
  bz.cpp
  debug.cpp
  lattice.cpp
  neighbours.cpp
  pointgroup.cpp
  spg_database.cpp
  symmetry.cpp
  trellis.cpp
)
foreach(CXX_TARGET ${CXX_TARGETS})
  target_sources(${CXX_TARGET} PRIVATE ${CXX_SOURCES})
  target_include_directories(${CXX_TARGET} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
endforeach(CXX_TARGET)
