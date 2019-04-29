add_custom_target (
        write_version_info
        COMMAND python ${CMAKE_SOURCE_DIR}/write_version_info.py ${CMAKE_BINARY_DIR}/version_info.h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

include_directories(${CMAKE_BINARY_DIR})

execute_process(COMMAND python -c "from write_version_info import print_version_number; print_version_number()"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE SYMBZ_VERSION)

# message(${SYMBZ_VERSION})
