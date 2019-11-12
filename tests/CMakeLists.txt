file(GLOB CATCH_SOURCES *_test.cpp)
target_sources(catch2 PRIVATE ${CATCH_SOURCES})
target_include_directories(catch2 PUBLIC ${CMAKE_CURRENT_LIST_DIR})
