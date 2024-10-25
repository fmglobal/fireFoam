cmake_minimum_required(VERSION 3.12)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(fireFoam LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-place build not allowed! Please use a separate build directory.")
endif ()

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS
        "No build type selected, default to Release. Use -DCMAKE_BUILD_TYPE=Debug for debug mode")
    set(CMAKE_BUILD_TYPE "Release")
endif()
message("-- Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_INCLUDE_CURRENT_DIR OFF)
#set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_RULE_MESSAGES OFF)

add_compile_definitions(CMAKE_VERSION_PRINT=TRUE)
set(FF_PKG_DIR ${CMAKE_SOURCE_DIR}/fireFoam/packages)
set(FF_LIB_DIR ${CMAKE_SOURCE_DIR}/fireFoam/lib)
set(FF_SRC_DIR ${CMAKE_SOURCE_DIR}/fireFoam/solver)

find_package(Git)
if(GIT_FOUND)
    string(TIMESTAMP _time_stamp)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE _build_version
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    message(STATUS "GIT hash: ${_build_version}")
    message(STATUS "time stamp: ${_time_stamp}")
else()
    # TODO: fallback behaviors to be implemented
    message(FATAL_ERROR "GIT not found")
endif()

configure_file(${CMAKE_SOURCE_DIR}/cmake/git_version.H.in ${FF_SRC_DIR}/git_version.H @ONLY)

include(OpenFOAM)

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES "$ENV{FOAM_USER_APPBIN}")
set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES "$ENV{FOAM_USER_LIBBIN}")

set(CMAKE_INSTALL_RPATH "$ENV{FOAM_USER_LIBBIN}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
message(STATUS "fireFOAM lib path: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(STATUS "fireFOAM bin path: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "fireFOAM install runtime directory: $ENV{FOAM_USER_APPBIN}")
message(STATUS "fireFOAM install library directory: $ENV{FOAM_USER_LIBBIN}")


add_subdirectory(fireFoam)
