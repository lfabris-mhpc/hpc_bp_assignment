cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.0 CACHE STRING "Floating point number for the current version")
project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

add_library(ljmd_rt SHARED
    src/engine.c
    src/input_output.c
    src/utilities.c
    src/defs.c
)
target_include_directories(ljmd_rt PRIVATE include)
target_link_libraries(ljmd_rt PRIVATE m)

add_executable(ljmd.x main.c)
target_include_directories(ljmd.x PRIVATE include)
target_compile_options(ljmd.x PUBLIC -Wall -g -O3 -march=native -ffast-math -fomit-frame-pointer -fexpensive-optimizations)
target_compile_definitions(ljmd.x PRIVATE -DLJMD_VERSION=${LJMD_VERSION})

target_link_libraries(ljmd.x PRIVATE m)
target_link_libraries(ljmd.x PRIVATE ljmd_rt)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
	
	set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
	enable_language(CXX)
	include(GTest)
	
    add_subdirectory(tests)
endif()
