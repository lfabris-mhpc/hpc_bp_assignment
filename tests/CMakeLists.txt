cmake_minimum_required(VERSION 3.10)

find_package(MPI REQUIRED)

#regression tests
file(COPY ../examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../reference DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY regression.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME Regression-test-argon_108 COMMAND regression.sh "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ../ljmd.x ${MPIEXEC_POSTFLAGS}" argon_108)
add_test(NAME Regression-test-argon_2916 COMMAND regression.sh "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS}../ljmd.x ${MPIEXEC_POSTFLAGS}" argon_2916)

#unit tests
#utilities
add_executable(test_utilities.x test_utilities.cpp)
target_include_directories(test_utilities.x PRIVATE ../include)
target_link_libraries(test_utilities.x PRIVATE GTest::GTestMain GTest::GTest ljmd_rt)

add_test(NAME Unit-test-utilities COMMAND test_utilities.x)

#engine
add_executable(test_engine.x test_engine.cpp)
target_include_directories(test_engine.x PRIVATE ../include)
target_link_libraries(test_engine.x PRIVATE GTest::GTestMain GTest::GTest ljmd_rt MPI::MPI_CXX MPI::MPI_C)

add_test(NAME Unit-test-engine COMMAND echo "argon_2.inp" > test_engine.x)

