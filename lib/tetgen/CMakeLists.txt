if (MSVC)
  add_compile_options(/O2)
else()
  add_compile_options(-O3)
endif()

set(TETGEN_TARGETS _symbz test_cpp test_catch)

set(TETGEN_SOURCES tetgen.cxx predicates.cxx)
foreach(TETGEN_TARGET ${TETGEN_TARGETS})
  target_sources(${TETGEN_TARGET} PRIVATE ${TETGEN_SOURCES})
  target_include_directories(${TETGEN_TARGET} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
endforeach(TETGEN_TARGET)
