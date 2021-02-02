cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.0 CACHE STRING "Floating point number for the current version")
project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

find_package(OpenMP)
if(OpenMP_C_FOUND)
#target_link_libraries(ljmd.omp.x PUBLIC OpenMP::OpenMP_C)
endif()

find_package(MPI)
if(MPI_FOUND)
	message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
	message(STATUS "Libs: ${MPI_INCLUDE_DIRS} ${MPI_C} ${MPI_C}")

	add_library(ljmd_rt.mpi SHARED
			src/engine.c
			src/input_output.c
			src/utilities.c
			)
	target_compile_options(ljmd_rt.mpi PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer -D__MPI)
	target_include_directories(ljmd_rt.mpi PRIVATE include ${MPI_INCLUDE_PATH})
	target_link_libraries(ljmd_rt.mpi PRIVATE m MPI::MPI_C)

	add_executable(ljmd.mpi.x main.c)
	target_include_directories(ljmd.mpi.x PRIVATE include)
	target_compile_options(ljmd.mpi.x PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer -D__MPI)
	target_compile_definitions(ljmd.mpi.x PRIVATE -DLJMD_VERSION=${LJMD_VERSION})

	target_link_libraries(ljmd.mpi.x PRIVATE m ljmd_rt.mpi MPI::MPI_C)
endif()

add_library(ljmd_rt SHARED
    src/engine.c
    src/input_output.c
    src/utilities.c
)
target_include_directories(ljmd_rt PRIVATE include)
target_link_libraries(ljmd_rt PRIVATE m)

add_executable(ljmd.x main.c)
target_include_directories(ljmd.x PRIVATE include)
target_compile_options(ljmd.x PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer)
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
