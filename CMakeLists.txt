cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(VULKAN_SDK_VERSION "Version of Vulkan SDK to use" "1.4.321.0")

# Utility function for fetching dependencies in one line.
function(fetch_dependency AGIT_REPO AGIT_TAG ANAME)
    # Remaining args are variable=value pairs for cache options
    set(options)
    set(oneValueArgs)
    set(multiValueArgs)
    cmake_parse_arguments(FETCHDEP "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Declare the dependency
    FetchContent_Declare(
        ${ANAME}
        GIT_REPOSITORY ${AGIT_REPO}
        GIT_TAG        ${AGIT_TAG}
        SOURCE_DIR     ${CMAKE_SOURCE_DIR}/dep/${ANAME}-src
    )

    # Apply cache variables
    foreach(arg IN LISTS ARGN)
        string(REPLACE "=" ";" kv "${arg}")
        list(GET kv 0 key)
        list(GET kv 1 val)
        set(${key} ${val} CACHE BOOL "")
    endforeach()

    # Actually fetch it
    FetchContent_MakeAvailable(${ANAME})
endfunction()

include(FetchContent)

set(PROJECT_NAME "geodesy-math")

file(GLOB_RECURSE INC
    "inc/*.h"
)

file(GLOB_RECURSE SRC
    "src/*.h"
    "src/*.cpp"
    "src/*.c"
)

project(${PROJECT_NAME})

add_library(${PROJECT_NAME} ${INC} ${SRC})

# Add Compiler flags.
if (MSVC)
    # MSVC projects are multi-config, so use generator expressions
    # MSVC: /O2 = optimize, /fp:fast = relaxed floating point,
    # /arch:AVX2 = use AVX2 if available
    target_compile_options(${GEODESY_LIBRARY} PUBLIC 
        $<$<CONFIG:Debug,RelWithDebInfo>:>
        $<$<CONFIG:Release,MinSizeRel>:
            /O2 
            /fp:fast 
            /arch:AVX2
        >
    )
    message(STATUS "Geodesy: MSVC SIMD flags enabled")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # GCC or Clang: optimize, relaxed FP, enable AVX2/FMA
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        # target_compile_options(${GEODESY_LIBRARY} PUBLIC -O0 -g)
    else()
        target_compile_options(${GEODESY_LIBRARY} PUBLIC 
            -O3 
            -ffast-math 
            -mavx2 
            -mfma
        )
    endif()
    message(STATUS "Geodesy: GCC/Clang SIMD flags enabled")
else()
    message(WARNING "Geodesy: No SIMD flags set for this compiler")
endif()

# Set Include Directories
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc/)
