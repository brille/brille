find_package(OpenMP REQUIRED) # Change code to support missing OpenMP?
if(OpenMP_CXX_FOUND)
    foreach(OMP_TARGET IN LISTS CXX_TARGETS)
        # message(STATUS "Adding OpenMP library to ${OMP_TARGET}")
        target_link_libraries(${OMP_TARGET} PUBLIC OpenMP::OpenMP_CXX)
    endforeach(OMP_TARGET)
    if (MSVC AND MSVC_VERSION GREATER 1919)
        add_compile_options(/openmp:experimental) # this doesn't work
    endif()
else()
    # The macOS llvm-openmp does sets OpenMP_FOUND and defines OpenMP::OpenMP as the library target...
    foreach(OMP_TARGET IN LISTS CXX_TARGETS)
        target_link_libraries(${OMP_TARGET} PUBLIC OpenMP::OpenMP)
    endforeach()
endif()