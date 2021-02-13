cmake_minimum_required(VERSION 3.10)

find_package(MPI REQUIRED)

#regression tests
file(COPY ../examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../reference DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY regression.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(mpi_cmd_serial "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} -x OMP_NUM_THREADS=1 ../ljmd.x ${MPIEXEC_POSTFLAGS}")
set(mpi_cmd_multip_omp "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} --map-by node:PE=${MPIEXEC_MAX_NUMPROCS} -x OMP_NUM_THREADS=${MPIEXEC_MAX_NUMPROCS} ../ljmd.x ${MPIEXEC_POSTFLAGS}")
set(mpi_cmd_multip_mpi "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} -x OMP_NUM_THREADS=1 ../ljmd.x ${MPIEXEC_POSTFLAGS}")
set(mpi_cmd_multip_hybrid "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} --map-by node:PE=2 -x OMP_NUM_THREADS=2 ../ljmd.x ${MPIEXEC_POSTFLAGS}")

add_test(NAME Regression-test-serial-argon_108 COMMAND regression.sh "${mpi_cmd_serial}" argon_108)
add_test(NAME Regression-test-multip-omp-argon_108 COMMAND regression.sh "${mpi_cmd_multip_omp}" argon_108)
add_test(NAME Regression-test-multip-mpi-argon_108 COMMAND regression.sh "${mpi_cmd_multip_mpi}" argon_108)
add_test(NAME Regression-test-multip-hybrid-argon_108 COMMAND regression.sh "${mpi_cmd_multip_hybrid}" argon_108)

#should add option for executing the larger regression tests
#add_test(NAME Regression-test-serial-argon_2916 COMMAND regression.sh "${mpi_cmd_serial}" argon_2916)
#add_test(NAME Regression-test-multip-omp-argon_2916 COMMAND regression.sh "${mpi_cmd_multip_omp}" argon_2916)
#add_test(NAME Regression-test-multip-mpi-argon_2916 COMMAND regression.sh "${mpi_cmd_multip_mpi}" argon_2916)
#add_test(NAME Regression-test-multip-hybrid-argon_2916 COMMAND regression.sh "${mpi_cmd_multip_hybrid}" argon_2916)

#unit tests
#utilities
add_executable(test_utilities.x test_utilities.cpp)
target_include_directories(test_utilities.x PRIVATE ../include)
target_link_libraries(test_utilities.x PRIVATE GTest::GTestMain GTest::GTest ljmd_rt)

add_test(NAME Unit-test-utilities COMMAND test_utilities.x)

#engine
add_executable(test_engine.x test_engine.cpp)
target_include_directories(test_engine.x PRIVATE ../include)
target_link_libraries(test_engine.x PRIVATE GTest::GTestMain GTest::GTest ljmd_rt)

add_test(NAME Unit-test-engine COMMAND echo "argon_2.inp" > test_engine.x)

#defs and input_output
add_executable(test_mdsys_alloc.x test_mdsys_alloc.cpp)
target_include_directories(test_mdsys_alloc.x PRIVATE ../include)
target_link_libraries(test_mdsys_alloc.x PRIVATE GTest::GTestMain GTest::GTest ljmd_rt)

add_test(NAME Unit-test-mdsys-alloc COMMAND echo "argon_2.inp" > test_mdsys_alloc.x)

#engine MPI functions
add_executable(test_engine_mpi.x test_engine_mpi.cpp)
target_include_directories(test_engine_mpi.x PRIVATE ../include)
target_link_libraries(test_engine_mpi.x PRIVATE GTest::GTestMain GTest::GTest ljmd_rt MPI::MPI_CXX MPI::MPI_C)

add_test(NAME Unit-test-engine-serial COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_engine_mpi.x)
add_test(NAME Unit-test-engine-multip-omp COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 --map-by node:PE=2 -x OMP_NUM_THREADS=${MPIEXEC_MAX_NUMPROCS} test_engine_mpi.x)
add_test(NAME Unit-test-engine-multip-mpi COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} --map-by node:PE=1 -x OMP_NUM_THREADS=1 test_engine_mpi.x)
add_test(NAME Unit-test-engine-multip-hybrid COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 --map-by node:PE=2 -x OMP_NUM_THREADS=2 test_engine_mpi.x)
