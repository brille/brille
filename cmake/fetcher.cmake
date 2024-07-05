include(FetchContent)

function(git_fetch package version source)
    FetchContent_Declare(${package} GIT_REPOSITORY ${source} GIT_TAG v${version})
    FetchContent_MakeAvailable(${package})
    set(${package}_FOUND ON PARENT_SCOPE)
    set("${package}_SOURCE_DIR" "${${package}_SOURCE_DIR}" PARENT_SCOPE)
    set("${package}_BINARY_DIR" "${${package}_BINARY_DIR}" PARENT_SCOPE)
endfunction()
