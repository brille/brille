# if(CMAKE_CXX_CPPCHECK)
#   list(APPEND CMAKE_CXX_CPPCHECK "--suppress=*:${CMAKE_CURRENT_LIST_DIR}/*")
# endif()

if(BRILLE_BUILD_TESTING)
  target_include_directories(catch2 PUBLIC ${CMAKE_CURRENT_LIST_DIR})
endif()

add_subdirectory(tetgen)
