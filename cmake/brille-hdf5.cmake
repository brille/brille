# HighFive is used to simplify the interface with HDF5:
# Override some default settings
if (BRILLE_HDF5)
    set(HIGHFIVE_USE_BOOST OFF)
    set(HIGHFIVE_UNIT_TESTS OFF)
    set(HIGHFIVE_EXAMPLES OFF)
    set(HIGHFIVE_BUILD_DOCS OFF)
    git_fetch(highfive ${MINIMUM_HIGHFIVE_VERSION} ${FETCH_HIGHFIVE_REPO} ${REQUIRE_SYSTEM_HIGHFIVE})
    foreach(HF_TARGET IN LISTS CXX_TARGETS)
        # message(STATUS "Adding HighFive library to ${HF_TARGET}")
        target_link_libraries(${HF_TARGET} PUBLIC HighFive)
    endforeach()
endif()