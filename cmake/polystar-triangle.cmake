
include(FetchContent)
FetchContent_Declare(
        Triangle
        GIT_REPOSITORY https://github.com/wo80/Triangle.git
        GIT_TAG        d3d0ccc94789e7e760f71de568d9d605127bb954
        SOURCE_SUBDIR  src
)
FetchContent_MakeAvailable(Triangle)
add_compile_definitions(Triangle PUBLIC NO_FILE_IO CDT_ONLY)

foreach(TRI_TARGET IN LISTS CXX_TARGETS)
    # message(STATUS "Adding HighFive library to ${HF_TARGET}")
    target_link_libraries(${TRI_TARGET} PUBLIC Triangle::triangle-api)
endforeach()