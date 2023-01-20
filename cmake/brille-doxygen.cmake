# first we can indicate the documentation build as an option (default OFF)
option(RUN_DOXYGEN_WITH_ALL "Build Doxygen documentation without specifying 'docs' target" OFF)
option(USE_DOXYGEN "Look for and use Doxygen to build documentation" OFF)

# check if Doxygen is installed
if (USE_DOXYGEN)
    find_package(Doxygen QUIET)
    if (DOXYGEN_FOUND)
        # set input and output files
        set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        if(BUILD_DOC)
            add_custom_target( docs ALL
                    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating documentation with Doxygen"
                    VERBATIM )
        else()
            add_custom_target( docs
                    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generate documentation using target 'docs'"
                    VERBATIM )
        endif()
    else (DOXYGEN_FOUND)
        message(STATUS "Install Doxygen to build documentation")
    endif (DOXYGEN_FOUND)
endif()