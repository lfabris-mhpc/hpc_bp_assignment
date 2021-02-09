cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.0 CACHE STRING "Floating point number for the current version")

project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

#src files
list(APPEND src/engine.c src/defs.c src/input_output.c src/utilities.c)

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

option(BUILD_CompilerOpt "Build with compiler optimization" ON)

option(BUILD_AlgorithmicOpt "Build with algorithmic optimization" ON)

find_package(OpenMP)
option(BUILD_OpenMP "Build OpenMP" ${OpenMP_C_FOUND})

find_package(MPI)
option(BUILD_MPI "Build MPI" ${MPI_FOUND})

if(MPI_FOUND AND OpenMP_C_FOUND)
	option(BUILD_HYBRID "Build hybrid" ON)
endif()

#actual targets
message(STATUS "Setting basic serial targets")
add_library(ljmd_rt SHARED 
	src/engine.c
	src/input_output.c
	src/utilities
	src/defs.c
)
target_include_directories(ljmd_rt PRIVATE include)
target_link_libraries(ljmd_rt PRIVATE m)

add_executable(ljmd.x main.c)
target_include_directories(ljmd.x PRIVATE include)
target_compile_options(ljmd.x PRIVATE ${ccflags})

target_link_libraries(ljmd.x PRIVATE m ljmd_rt)

if(BUILD_CompilerOpt)
	add_library(ljmd_rt.copt SHARED 
		src/engine.c
		src/input_output.c
		src/utilities.c
		src/defs.c
		)
	#list(APPEND lccflags ${ccflags} ${ccflags_copt})
	
	target_include_directories(ljmd_rt.copt PRIVATE include)
	target_link_libraries(ljmd_rt.copt PRIVATE m)

	add_executable(ljmd.copt.x main.c)
	target_include_directories(ljmd.copt.x PRIVATE include)
	target_compile_options(ljmd.copt.x PRIVATE ${ccflags} ${ccflags_copt})

	target_link_libraries(ljmd.copt.x PRIVATE m ljmd_rt.copt)
endif()

if(BUILD_AlgorithmicOpt)
	add_library(ljmd_rt.aopt SHARED 
		src/engine.c
		src/input_output.c
		src/utilities.c
		src/defs.c
		)
	#list(APPEND lccflags ${ccflags} ${ccflags_copt} ${ccflags_aopt})
	
	target_include_directories(ljmd_rt.aopt PRIVATE include)
	target_link_libraries(ljmd_rt.aopt PRIVATE m)

	add_executable(ljmd.aopt.x main.c)
	target_include_directories(ljmd.aopt.x PRIVATE include)
	target_compile_options(ljmd.aopt.x PRIVATE ${ccflags} ${ccflags_copt} ${ccflags_aopt})

	target_link_libraries(ljmd.aopt.x PRIVATE m ljmd_rt.aopt)
endif()

if(BUILD_OpenMP)
	message(STATUS "Setting up OpenMP targets")
	
	add_library(ljmd_rt.omp SHARED 
		src/engine.c
		src/input_ouput.c
		src/utilities.c
		src/defs.c
		)
	#list(APPEND lccflags ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_omp})
	
	target_compile_options(ljmd_rt.omp PRIVATE ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_omp})
	target_include_directories(ljmd_rt.omp PRIVATE include)
	target_link_libraries(ljmd_rt.omp PRIVATE m OpenMP::OpenMP_C)

	add_executable(ljmd.omp.x main.c)
	target_include_directories(ljmd.omp.x PRIVATE include)
	target_compile_options(ljmd.omp.x PRIVATE ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_omp})

	target_link_libraries(ljmd.omp.x PRIVATE m ljmd_rt.omp OpenMP::OpenMP_C)
endif()

if(BUILD_MPI)
	message(STATUS "Setting up MPI targets")
	#message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
	#message(STATUS "Libs: ${MPI_INCLUDE_DIRS} ${MPI_C}")
	
	add_library(ljmd_rt.mpi SHARED 
		src/engine.c
		src/input_ouput.c
		src/utilities.c
		src/defs.c
		)
	#list(APPEND lccflags ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_mpi})
	
	target_compile_options(ljmd_rt.mpi PRIVATE ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_mpi})
	target_include_directories(ljmd_rt.mpi PRIVATE include ${MPI_INCLUDE_PATH})
	target_link_libraries(ljmd_rt.mpi PRIVATE m MPI::MPI_C)

	add_executable(ljmd.mpi.x main.c)
	target_include_directories(ljmd.mpi.x PRIVATE include)
	target_compile_options(ljmd.mpi.x PRIVATE ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_mpi})

	target_link_libraries(ljmd.mpi.x PRIVATE m ljmd_rt.mpi MPI::MPI_C)
endif()

if(BUILD_HYBRID)
	message(STATUS "Setting up hybrid MPI+OpenMP targets")
	
	add_library(ljmd_rt.hybrid SHARED 
		src/engine.c
		src/input_output.c
		src/utilities.c
		src/defs.c
		)
	#list(APPEND lccflags ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_mpi} ${ccflags_omp})

	target_compile_options(ljmd_rt.hybrid PRIVATE ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_mpi} ${ccflags_omp})
	target_include_directories(ljmd_rt.hybrid PRIVATE include ${MPI_INCLUDE_PATH})
	target_link_libraries(ljmd_rt.hybrid PRIVATE m MPI::MPI_C OpenMP::OpenMP_C)

	add_executable(ljmd.hybrid.x main.c)
	target_include_directories(ljmd.hybrid.x PRIVATE include)
	target_compile_options(ljmd.hybrid.x PRIVATE ${ccflags} ${ccflags_copt} ${ccflags_aopt} ${ccflags_mpi} ${ccflags_omp})

	target_link_libraries(ljmd.hybrid.x PRIVATE m ljmd_rt.hybrid MPI::MPI_C OpenMP::OpenMP_C)
endif()

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
	
	set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
	enable_language(CXX)
	include(GTest)
	
    add_subdirectory(tests)
endif()
