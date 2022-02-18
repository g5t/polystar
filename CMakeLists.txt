cmake_minimum_required(VERSION 3.18.2) 
cmake_policy(SET CMP0076 NEW) # Ensure target_sources converts relative paths

project(polystar)
set(BRILLE_LIBRARY_TARGET polystar)
set(BRILLE_PYTHON_MODULE _polystar)
set(BRILLE_SINGLE_HEADER polystar.h) # must match template file in project root
SET(BRILLE_LIB_DESTINATION lib)
SET(BRILLE_BIN_DESTINATION bin)
SET(BRILLE_INCLUDE_DESTINATION include)
set(TESTING_TARGET tester)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Make the availability of testing optional
option(BRILLE_BUILD_TESTING "Build tests for polystar" ON)
# Allow system pybind11 to be required
option(REQUIRE_SYSTEM_PYBIND11 "Never attempt to fetch pybind11" OFF)
mark_as_advanced(REQUIRE_SYSTEM_PYBIND11)
# Allow system catch2 to be required
option(REQUIRE_SYSTEM_CATCH2 "Never attempt to fetch catch2" OFF)
mark_as_advanced(REQUIRE_SYSTEM_CATCH2)
# Allow the log level to be set by cmake
set(BRILLE_LOGLEVEL "INFO" CACHE STRING "Emit log messages to standard out (DEBUG|VERBOSE)")
# Special option for profiling runs
option(BRILLE_PROFILING "Emit profiling output to standard out" OFF)

# Define the minimum version of libraries needed
set(MINIMUM_PYBIND11_VERSION 2.8.1)
set(FETCH_PYBIND11_REPO https://github.com/pybind/pybind11)
set(MINIMUM_CATCH2_VERSION 2.11.3)
set(FETCH_CATCH2_REPO https://github.com/catchorg/Catch2)

# Optional support of HDF5-based IO
option(BRILLE_HDF5 "Add HDF5 file-based IO" OFF)
option(REQUIRE_SYSTEM_HIGHFIVE "Never attempt to fetch HighFive" OFF)
mark_as_advanced(REQUIRE_SYSTEM_HIGHFIVE)
set(MINIMUM_HIGHFIVE_VERSION 2.3.1)
SET(FETCH_HIGHFIVE_REPO https://github.com/BlueBrain/HighFive)

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
else()
    # lots of warnings -- add -Werror for  all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CMAKE_MACOSX_RPATH 1)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

if (${BRILLE_LOGLEVEL} STREQUAL "VERBOSE")
  message(STATUS "Verbose logging emitted at runtime")
  add_definitions(-DVERBOSE)
else()
  if (${BRILLE_LOGLEVEL} STREQUAL "DEBUG")
    message(STATUS "Debug logging emitted at runtime")
    add_definitions(-DDEBUG)
  else()
    message(STATUS "Informational logging emitted at runtime")
  endif()
endif()
if (BRILLE_PROFILING)
  message(STATUS "Profiling output emitted at runtime")
  add_definitions(-DPROFILING)
endif (BRILLE_PROFILING)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (PYTHON_EXECUTABLE)
  # Ensure the provided Python interpreter is used
  set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
endif()
# With modern CMake, this find_package forces pybind11 to use FindPython instead of its custom tools
find_package(Python3 COMPONENTS Interpreter Development)

include(fetcher.cmake)

# Attempt to find catch2 to handle C++ testing
if(BRILLE_BUILD_TESTING)
  git_fetch(catch2 ${MINIMUM_CATCH2_VERSION} ${FETCH_CATCH2_REPO} ${REQUIRE_SYSTEM_CATCH2})
  list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/contrib")
else()
  # Since no testing is to be built, fake that we've found Catch2.
  set(Catch2_FOUND ON)
endif()
# Attempt to find pybind11 to handle CPython bindings
git_fetch(pybind11 ${MINIMUM_PYBIND11_VERSION} ${FETCH_PYBIND11_REPO} ${REQUIRE_SYSTEM_PYBIND11})

# HighFive is used to simplify the interface with HDF5:
# Override some default settings
if (BRILLE_HDF5)
  set(HIGHFIVE_USE_BOOST OFF)
  set(HIGHFIVE_UNIT_TESTS OFF)
  set(HIGHFIVE_EXAMPLES OFF)
  set(HIGHFIVE_BUILD_DOCS OFF)
  git_fetch(highfive ${MINIMUM_HIGHFIVE_VERSION} ${FETCH_HIGHFIVE_REPO} ${REQUIRE_SYSTEM_HIGHFIVE})
  if (highfive_FOUND)
		message(STATUS "Adding HighFive precompiler macro flag")
    add_definitions(-DUSE_HIGHFIVE)
  endif()
endif(BRILLE_HDF5)

# Read the version of polystar
execute_process(COMMAND ${Python3_EXECUTABLE} -c "from version_info import version_number; version_number()"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BRILLE_VERSION)
# So that we can print it to the console along with the specified build type
message(STATUS "Build polystar v${BRILLE_VERSION} with type ${CMAKE_BUILD_TYPE}")

# Create a header file with information about the last time tye python module was built
# This file will be updated even if no other source files are modified but it does not
# cause relinking the module (if it is the *only* file which changes)
add_custom_target(write_version_info
  COMMAND ${Python3_EXECUTABLE} version_info.py src/version.hpp
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# Create a single header by contatenating all headers in src/
add_custom_target(single_header
  COMMAND ${Python3_EXECUTABLE} acme.py ${BRILLE_SINGLE_HEADER} -o ${CMAKE_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  DEPENDS write_version_info
)

add_custom_target(library DEPENDS single_header ${BRILLE_LIBRARY_TARGET})


# We will always build the python module
list(APPEND CXX_TARGETS ${BRILLE_PYTHON_MODULE})

if(BRILLE_BUILD_TESTING)
  enable_testing() # allows registration of Python tests in wrap/
endif()

# Target for python module
pybind11_add_module(${BRILLE_PYTHON_MODULE} MODULE)
add_subdirectory(wrap)

if(BRILLE_BUILD_TESTING)
  list(APPEND CXX_TARGETS ${TESTING_TARGET}) # Include the C++ test target
  # target for Catch2 based tests
  add_executable(${TESTING_TARGET})
endif()
add_subdirectory(lib)

# target for C++ shared library
add_library(${BRILLE_LIBRARY_TARGET} src/hdf_interface.hpp)
list(APPEND CXX_TARGETS ${BRILLE_LIBRARY_TARGET})

# add the dependencies and include directories for all CXX targets:
add_subdirectory(src)  # important

if(BRILLE_BUILD_TESTING)
  target_link_libraries(${TESTING_TARGET} PUBLIC Catch2::Catch2)
  include(CTest)
  include(Catch)
  catch_discover_tests(${TESTING_TARGET})
endif()

# OpenMP support:
find_package(OpenMP REQUIRED) # Change code to support missing OpenMP?
if(OpenMP_CXX_FOUND)
  foreach(OMP_TARGET IN LISTS CXX_TARGETS)
    target_link_libraries(${OMP_TARGET} PUBLIC OpenMP::OpenMP_CXX)
  endforeach(OMP_TARGET)
  if (MSVC AND MSVC_VERSION GREATER 1919)
    add_compile_options(/openmp:experimental) # this doesn't work
  endif()
endif()

if(BRILLE_HDF5 AND highfive_FOUND)
  foreach(HF_TARGET IN LISTS CXX_TARGETS)
    message(STATUS "Adding HighFive library to ${HF_TARGET}")
    target_link_libraries(${HF_TARGET} PUBLIC HighFive)
  endforeach()
endif()


# if(MSVC)
# 	if(CMAKE_BUILD_TYPE STREQUAL "Matlab")
# 		find_package(Matlab)
# 		target_link_libraries(_polystar PUBLIC "${Matlab_ROOT_DIR}/bin/win64/libiomp5md.lib")
# 		set_property(TARGET _polystar APPEND PROPERTY LINK_FLAGS /nodefaultlib:vcomp)
# 	endif()
# else()
# endif()

# first we can indicate the documentation build as an option (default OFF)
option(BUILD_DOC "Build documentation" OFF)
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

install(
  TARGETS ${BRILLE_LIBRARY_TARGET}
  ARCHIVE DESTINATION ${BRILLE_LIB_DESTINATION}
  LIBRARY DESTINATION ${BRILLE_LIB_DESTINATION}
  RUNTIME DESTINATION ${BRILLE_BIN_DESTINATION}
)
install(
  FILES "${CMAKE_BINARY_DIR}/${BRILLE_SINGLE_HEADER}"
  DESTINATION ${BRILLE_INCLUDE_DESTINATION}
)