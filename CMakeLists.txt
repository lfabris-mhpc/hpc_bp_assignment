cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.1 CACHE STRING "Floating point number for the current version")

project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

#src files
list(APPEND src_files src/engine.c src/defs.c src/input_output.c src/utilities.c)

list(APPEND ccflags "-Wall" "-Wextra" "-DLJMD_VERSION=${LJMD_VERSION}")
list(APPEND ccflags_copt "-g" "-O3" "-march=native" "-DNDEBUG" "-ffast-math" "-fomit-frame-pointer" "-fexpensive-optimizations")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_library(ljmd_rt SHARED ${src_files})
target_include_directories(ljmd_rt PRIVATE include)
target_link_libraries(ljmd_rt PRIVATE m MPI::MPI_C)

add_executable(ljmd.x main.c)
target_include_directories(ljmd.x PRIVATE include)
target_compile_options(ljmd.x PRIVATE ${ccflags} ${ccflags_opt})

target_link_libraries(ljmd.x PRIVATE m ljmd_rt MPI::MPI_C)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
	
	set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
	enable_language(CXX)
	include(GTest)
	
    add_subdirectory(tests)
endif()
