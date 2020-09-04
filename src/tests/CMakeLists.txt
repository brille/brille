file(GLOB CATCH_SOURCES *_test.cpp)
target_sources(${TESTING_TARGET} PRIVATE ${CATCH_SOURCES})
target_include_directories(${TESTING_TARGET} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
