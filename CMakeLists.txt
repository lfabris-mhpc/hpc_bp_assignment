cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.0 CACHE STRING "Floating point number for the current version")

project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

#src files
list(APPEND src_files src/engine.c src/defs.c src/input_output.c src/utilities.c)

#basic compiler flags
list(APPEND ccflags "-Wall" "-Wextra" "-DLJMD_VERSION=${LJMD_VERSION}")
#compiler-optimization flags
list(APPEND ccflags_copt "-g" "-O3" "-march=native" "-DNDEBUG" "-ffast-math" "-fomit-frame-pointer" "-fexpensive-optimizations")
#algorithmic-optimization enabling flags
list(APPEND ccflags_aopt "-D__OPT")
#OpenMP-enabling flags
list(APPEND ccflags_omp "-fopenmp")
#MPI-enabling flags
list(APPEND ccflags_mpi "-D__MPI")

find_package(OpenMP)

#actual targets
message(STATUS "Setting basic serial targets")
add_library(ljmd_rt SHARED ${src_files})
target_include_directories(ljmd_rt PRIVATE include)
target_compile_options(ljmd_rt PRIVATE ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_omp})
target_link_libraries(ljmd_rt PRIVATE m OpenMP::OpenMP_C)

add_executable(ljmd.x main.c)
target_include_directories(ljmd.x PRIVATE include)
target_compile_options(ljmd_rt PRIVATE ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_omp})

target_link_libraries(ljmd.x PRIVATE m ljmd_rt OpenMP::OpenMP_C) 

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
	
	set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
	enable_language(CXX)
	include(GTest)
	
    add_subdirectory(tests)
endif()
