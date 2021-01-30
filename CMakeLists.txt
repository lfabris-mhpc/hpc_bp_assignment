cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.0 CACHE STRING "Floating point number for the current version")
project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

add_executable(ljmd.x
	main.c
)

add_library(ljmd_rt SHARED
    src/engine.c
    src/input_output.c
    src/utilities.c
)

target_include_directories(ljmd.x PRIVATE include)
target_include_directories(ljmd_rt PRIVATE include)

target_compile_options(ljmd.x PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer)
target_compile_definitions(ljmd.x PRIVATE -DLJMD_VERSION=${LJMD_VERSION})

target_link_libraries(ljmd.x PRIVATE m)
target_link_libraries(ljmd_rt PRIVATE m)
target_link_libraries(ljmd.x PRIVATE ljmd_rt)
