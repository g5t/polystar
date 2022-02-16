# Directory specific compilation flags
if (MSVC)
    # warning level 4 -- add /WX for all warnings as errors
    # /Zc__cplusplus ensures the preprocessor macro __cplusplus is set
    add_compile_options(/W4 /Zc__cplusplus)
else()
    # lots of warnings -- add -Werror for  all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic -Wmissing-braces)
endif()

# if(CMAKE_CXX_CPPCHECK)
#   list(APPEND CMAKE_CXX_CPPCHECK
#     "--suppress=compareBoolExpressionWithInt:${CMAKE_CURRENT_LIST_DIR}/arrayvector.hpp:349"
#   )
# endif()

message(STATUS "Using interpreter: ${Python3_EXECUTABLE}")

#set(CXX_SOURCES
list(APPEND CXX_SOURCES
  comparisons.cpp
  debug.cpp
)
if(BRILLE_HDF5 AND highfive_FOUND)
    message(STATUS "Adding hdf_interface.cpp")
    list(APPEND CXX_SOURCES hdf_interface.cpp)
endif()

foreach(CXX_TARGET IN LISTS CXX_TARGETS)
  target_sources(${CXX_TARGET} PRIVATE ${CXX_SOURCES})
  target_include_directories(${CXX_TARGET} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
endforeach()

if(BRILLE_BUILD_TESTING)
  add_subdirectory(tests)
endif()
