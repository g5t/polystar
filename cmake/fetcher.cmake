include(FetchContent)

function(git_fetch package version source required)
    if (${required})
        find_package(${package} ${version} REQUIRED)
    else()
        find_package(${package} ${version} QUIET)
    endif()

    if (${${package}_FOUND})
        message(STATUS "Found system ${package}")
    else()
        message(STATUS "Fetch ${package} ${version} from ${source}")
        FetchContent_Declare(${package} GIT_REPOSITORY ${source} GIT_TAG v${version})
        FetchContent_GetProperties(${package})
        if (NOT "${package}_POPULATED")
            FetchContent_Populate(${package})
            add_subdirectory("${${package}_SOURCE_DIR}" "${${package}_BINARY_DIR}")
        endif()
        set(${package}_FOUND ON PARENT_SCOPE)
        set("${package}_SOURCE_DIR" "${${package}_SOURCE_DIR}" PARENT_SCOPE)
        set("${package}_BINARY_DIR" "${${package}_BINARY_DIR}" PARENT_SCOPE)
    endif()
endfunction()

function(git_fetch_hash package source hash subdir)
    find_package(${package} QUIET)
    if (${${package}_FOUND})
        message(STATUS "Found ${package}")
    else()
        message(STATUS "Fetch ${package} ${hash} from ${source}")
        FetchContent_Declare(${package} GIT_REPOSITORY ${source} GIT_TAG ${hash})
        FetchContent_GetProperties(${package})
        if (NOT "${package}_POPULATED")
            FetchContent_Populate(${package})
            add_subdirectory("${${package}_SOURCE_DIR}" "${${package}_BINARY_DIR}")
        endif()
        set(${package}_FOUND ON PARENT_SCOPE)
        set("${package}_SOURCE_DIR" "${${package}_SOURCE_DIR}" PARENT_SCOPE)
        set("${package}_BINARY_DIR" "${${package}_BINARY_DIR}" PARENT_SCOPE)
    endif()
endfunction()