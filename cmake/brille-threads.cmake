set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
foreach(TARGET IN LISTS CXX_TARGETS)
    target_link_libraries(${TARGET} PRIVATE Threads::Threads)
    if (LINUX)
        # Set thread stack size to 8MB on Linux (the apparent default for glibc)
        # musllibc has default 128k https://wiki.musl-libc.org/functional-differences-from-glibc
        # which is a problem for using tetgen since it has many large stack arrays.
        target_link_options(${TARGET} PRIVATE "-Wl,-z,stack-size=8388608")
    endif()
endforeach()