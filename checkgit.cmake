set(CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR})

if (NOT DEFINED pre_configure_dir)
  set(pre_configure_dir ${CMAKE_CURRENT_LIST_DIR})
endif()

if (NOT DEFINED post_configure_dir)
  set(post_configure_dir ${CMAKE_BINARY_DIR}/gen)
endif()

if (NOT EXISTS ${post_configure_dir})
  message(STATUS "Create ${post_configure_dir}")
  file(MAKE_DIRECTORY ${post_configure_dir})
else()
  message(STATUS "Use ${post_configure_dir}")
endif()

set(pre_configure_file ${pre_configure_dir}/version.hpp.in)
set(post_configure_file ${post_configure_dir}/version.hpp)

function(checkGitWrite git_hash)
  file(WRITE ${CMAKE_BINARY_DIR}/git-state ${git_hash})
endfunction()

function(checkGitRead git_hash)
  if (EXISTS ${CMAKE_BINARY_DIR}/git-state)
    file(STRINGS ${CMAKE_BINARY_DIR}/git-state CONTENT)
    LIST(GET CONTENT 0 var)
    set(${git_hash} ${var} PARENT_SCOPE)
  endif()
endfunction()

function(checkGitVersion git_version)

  execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(
    COMMAND python -c "from datetime import datetime; print(datetime.now().isoformat(timespec='minutes'))"
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_CONFIGURE_TIME
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(
    COMMAND python -c "import setuptools_scm as s; print(s.get_version())"
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(
    COMMAND python -c "import setuptools_scm as s; print('.'.join(s.get_version().split('.')[:3]))"
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_SAFE_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(
    COMMAND python -c "import platform; print(platform.node())"
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_HOSTNAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  set(${git_version} ${GIT_VERSION} PARENT_SCOPE)
  checkGitRead(GIT_HASH_CACHE)

  if (NOT DEFINED GIT_HASH_CACHE)
    set(GIT_HASH_CACHE "INVALID")
  endif()

  if (NOT ${GIT_HASH} STREQUAL ${GIT_HASH_CACHE} OR NOT EXISTS ${post_configure_file})
    checkGitWrite(${GIT_HASH})
    configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
  endif()

endfunction()

function(checkGitSetup)
  add_custom_target(AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
    -DRUN_CHECK_GIT_VERSION=1
    -Dpre_configure_dir=${pre_configure_dir}
    -Dpost_configure_dir=${post_configure_dir}
    -DGIT_HASH_CACHE=${GIT_HASH_CACHE}
    -P ${CURRENT_LIST_DIR}/checkgit.cmake
    BYPRODUCTS ${post_configure_file}
  )
  checkGitVersion(GIT_VERSION)
  if (NOT DEFINED GIT_VERSION)
    set(GIT_VERSION "UNKNONW")
  endif()
  set(GIT_VERSION ${GIT_VERSION} PARENT_SCOPE)
  #add_library(git_version ${post_configure_dir}/version.hpp)
  #target_include_directories(git_version PUBLIC ${post_configure_dir})
  #add_dependencies(git_version AlwaysCheckGit)
  # checkGitVersion()
endfunction()

function(addGitVersion target)
  target_include_directories(${target} PUBLIC ${post_configure_dir})
  add_dependencies(${target} AlwaysCheckGit)
endfunction()


if (RUN_CHECK_GIT_VERSION)
  checkGitVersion(GIT_VERSION)
  if (NOT DEFINED GIT_VERSION)
    set(GIT_VERSION "UNKNOWN")
  endif()
endif()
