cmake_minimum_required(VERSION 3.10)

find_package(MPI REQUIRED)

#regression tests
file(COPY ../examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../reference DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY regression.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(mpi_cmd_serial "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ../ljmd.x ${MPIEXEC_POSTFLAGS}")
set(mpi_cmd_multip "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ../ljmd.x ${MPIEXEC_POSTFLAGS}")

add_test(NAME Regression-test-serial-argon_108 COMMAND regression.sh "${mpi_cmd_serial}" argon_108)
add_test(NAME Regression-test-multip-argon_108 COMMAND regression.sh "${mpi_cmd_multip}" argon_108)
#add_test(NAME Regression-test-serial-argon_2916 COMMAND regression.sh "${mpi_cmd_serial}" argon_2916)
#add_test(NAME Regression-test-multip-argon_2916 COMMAND regression.sh "${mpi_cmd}" argon_2916)

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

add_executable(test_engine_mpi.x test_engine_mpi.cpp)
target_include_directories(test_engine_mpi.x PRIVATE ../include)
target_link_libraries(test_engine_mpi.x PRIVATE GTest::GTestMain GTest::GTest ljmd_rt MPI::MPI_CXX MPI::MPI_C)

set(test_mpi_cmd_serial "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ../test_engine_mpi.x ${MPIEXEC_POSTFLAGS}")
set(test_mpi_cmd_multip "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ../test_engine_mpi.x ${MPIEXEC_POSTFLAGS}")

add_test(NAME Unit-test-engine-serial COMMAND mpirun --np 1 test_engine_mpi.x)
add_test(NAME Unit-test-engine-multip COMMAND mpirun --np 2 test_engine_mpi.x)
