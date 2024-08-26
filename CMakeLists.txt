cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0076 NEW) # Ensure target_sources converts relative paths
cmake_policy(SET CMP0017 NEW) # Prefer cmake's own files for include/find_package before CMAKE_MODULE_PATH

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_PROJECT_TOP_LEVEL_INCLUDES "${CMAKE_CURRENT_LIST_DIR}/cmake/conan_provider.cmake")

project(polystar)

set(POLYSTAR_LIBRARY_TARGET polystar)
set(POLYSTAR_PYTHON_MODULE _polystar)
set(POLYSTAR_PYTHON_DESTINATION polystar)
set(POLYSTAR_SINGLE_HEADER polystar.h) # must match template file in project root
SET(POLYSTAR_LIB_DESTINATION lib)
SET(POLYSTAR_BIN_DESTINATION bin)
SET(POLYSTAR_INCLUDE_DESTINATION include)
set(TESTING_TARGET tester)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make the availability of testing optional
option(POLYSTAR_BUILD_TESTING "Build tests for polystar" ON)

# Allow the log level to be set by cmake
set(POLYSTAR_LOGLEVEL "INFO" CACHE STRING "Emit log messages to standard out (DEBUG|VERBOSE)")
# Special option for profiling runs
option(POLYSTAR_PROFILING "Emit profiling output to standard out" OFF)

# Optional support of HDF5-based IO
# Include the HDF5 IO methods -- safe since configuring will fail if HDF5/HighFive is not available
add_definitions(-DUSE_HIGHFIVE)


# find_program(CMAKE_CXX_CPPCHECK NAMES cppcheck)
# if (CMAKE_CXX_CPPCHECK)
#     list(
#         APPEND CMAKE_CXX_CPPCHECK
#             "--enable=warning"
#             "--inconclusive"
#             "--force"
#             "--inline-suppr"
#             # "--template=gcc" # uncomment to get suppression error ids in brackets
#     )
# endif()

if (MSVC)
    # warning level 4 -- add /WX for all warnings as errors
    add_compile_options(/W4)
    # suppress MSVC warning C4996 about 'localtime' vs 'localtime_s'
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    # Allow UTF-8 identifiers https://stackoverflow.com/a/47704050
    add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
    add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")
else()
    # lots of warnings -- add -Werror for  all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CMAKE_MACOSX_RPATH OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    add_definitions(-DDEBUG)
endif()

if (${POLYSTAR_LOGLEVEL} STREQUAL "VERBOSE")
  message(STATUS "Verbose logging emitted at runtime")
  add_definitions(-DVERBOSE)
else()
  if (${POLYSTAR_LOGLEVEL} STREQUAL "DEBUG")
    message(STATUS "Debug logging emitted at runtime")
    add_definitions(-DDEBUG)
  else()
    message(STATUS "Informational logging emitted at runtime")
  endif()
endif()
if (POLYSTAR_PROFILING)
  message(STATUS "Profiling output emitted at runtime")
  add_definitions(-DPROFILING)
endif (POLYSTAR_PROFILING)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# With GCC 10+ the interproceedural optmization only adds to compilation time without improving performance
SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)

# Force pybind11 to use FindPython:
set(PYBIND11_FINDPYTHON ON)
# Ensure we can find/user the user-provided Python executable
if (PYTHON_EXECUTABLE)
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
endif()
find_package(Python 3 COMPONENTS Interpreter Development)
if (NOT DEFINED PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
endif()

find_package(Catch2 REQUIRED)
find_package(pybind11 REQUIRED)

# Read the version of polystar
include(checkgit)
checkGitSetup(POLYSTAR) # defines/creates version.hpp; and library 'git_version'
# So that we can print it to the console along with the specified build type
message(STATUS "Build polystar v${POLYSTAR_VERSION} with type ${CMAKE_BUILD_TYPE}")

if (NOT SKBUILD)
# Create a single header by contatenating all headers in src/
add_custom_target(single_header
  COMMAND ${Python_EXECUTABLE} acme.py ${POLYSTAR_SINGLE_HEADER} -o ${CMAKE_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  DEPENDS git_version
)
add_custom_target(library DEPENDS single_header ${POLYSTAR_LIBRARY_TARGET})
endif()

# We will always build the python module
list(APPEND CXX_TARGETS ${POLYSTAR_PYTHON_MODULE})

if(POLYSTAR_BUILD_TESTING)
  enable_testing() # allows registration of Python tests in wrap/
endif()

# Target for python module
pybind11_add_module(${POLYSTAR_PYTHON_MODULE} MODULE)
add_subdirectory(wrap)

if(POLYSTAR_BUILD_TESTING)
  list(APPEND CXX_TARGETS ${TESTING_TARGET}) # Include the C++ test target
  # target for Catch2 based tests
  add_executable(${TESTING_TARGET})
endif()

if(NOT SKBUILD)
# target for C++ shared library
add_library(${POLYSTAR_LIBRARY_TARGET} src/hdf_interface.hpp)
list(APPEND CXX_TARGETS ${POLYSTAR_LIBRARY_TARGET})
endif()

add_subdirectory(lib)
# add the dependencies and include directories for all CXX targets:
add_subdirectory(src)  # important

foreach(CXX_TARGET IN LISTS CXX_TARGETS)
  addGitVersion(${CXX_TARGET})
endforeach()

if(POLYSTAR_BUILD_TESTING)
  message(STATUS "Build testing target '${TESTING_TARGET}' and configure CTest")
  target_link_libraries(${TESTING_TARGET} PUBLIC Catch2::Catch2WithMain)
  include(CTest)
  include(Catch)
  catch_discover_tests(${TESTING_TARGET})
endif()

include(polystar-openmp)
include(polystar-hdf5)
include(polystar-triangle)
include(polystar-doxygen)

if (NOT SKBUILD)
install(
  TARGETS ${POLYSTAR_LIBRARY_TARGET}
  ARCHIVE DESTINATION ${POLYSTAR_LIB_DESTINATION}
  LIBRARY DESTINATION ${POLYSTAR_LIB_DESTINATION}
  RUNTIME DESTINATION ${POLYSTAR_BIN_DESTINATION}
)
install(
  FILES "${CMAKE_BINARY_DIR}/${POLYSTAR_SINGLE_HEADER}"
  DESTINATION ${POLYSTAR_INCLUDE_DESTINATION}
)
endif()

install(TARGETS ${POLYSTAR_PYTHON_MODULE} DESTINATION ${POLYSTAR_PYTHON_DESTINATION})
