file(GLOB CATCH_SOURCES *_test.cpp)
target_sources(test_brille PRIVATE ${CATCH_SOURCES})
target_include_directories(test_brille PUBLIC ${CMAKE_CURRENT_LIST_DIR})
