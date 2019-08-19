file(GLOB CATCH_SOURCES *_test.cpp)
target_sources(test_catch PRIVATE ${CATCH_SOURCES})
target_include_directories(test_catch PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(test_cpp PRIVATE test.cpp)
