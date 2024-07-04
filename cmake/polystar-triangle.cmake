
include(FetchContent)
FetchContent_Declare(
        Triangle
        GIT_REPOSITORY https://github.com/wo80/Triangle.git
        GIT_TAG        d3d0ccc94789e7e760f71de568d9d605127bb954
        SOURCE_SUBDIR  src
)
FetchContent_MakeAvailable(Triangle)

add_compile_definitions(Triangle PUBLIC NO_FILE_IO) # The API interface requires CDT_ONLY is not set to compile sources

foreach(TRI_TARGET IN LISTS CXX_TARGETS)
#    # Triangle::triangle-api uses shared libraries, which can be problematic for {delvewheel|delocate} to find
#    target_link_libraries(${TRI_TARGET} PRIVATE Triangle::triangle-api)

    # Triangle::triangle uses static libraries, which do not need to be found after building the module
    target_link_libraries(${TRI_TARGET} PUBLIC Triangle::triangle)
    # but _do_ need to be specified more-accurately since they show up in the build directory under Triangle
    target_include_directories(${TRI_TARGET}
	    PUBLIC "${Triangle_BINARY_DIR}/Triangle"
	    PRIVATE "${Triangle_SOURCE_DIR}/Triangle" "${Triangle_SOURCE_DIR}/Triangle/private"
    )
    # Since we _use_ the API interface, we must actually include the sources
    target_sources(${TRI_TARGET} PRIVATE
            "${Triangle_SOURCE_DIR}/Triangle/triangle_api.c"
            "${Triangle_SOURCE_DIR}/Triangle/triangle_helper.c"
    )
endforeach()
