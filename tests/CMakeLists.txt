cmake_minimum_required(VERSION 3.10)

#dry run
add_test(NAME Dry-run-test COMMAND echo "" > ../ljmd.x)

#regression tests
file(COPY ../examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../reference DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY regression.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME Regression-test COMMAND regression.sh)

#unit tests
#utilities
add_executable(test_utilities.x test_utilities.cpp)
target_include_directories(test_utilities.x PRIVATE ../include)
target_link_libraries(test_utilities.x PRIVATE GTest::GTestMain GTest::GTest ljmd_rt)

add_test(NAME Unit-test-utilities COMMAND test_utilities.x)

if(MPI_FOUND)
	file(COPY regression-MPI.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
	
	#add_test(NAME Regression-test-MPI COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ../ljmd.mpi.x}
	add_test(NAME Regression-test-MPI COMMAND regression-MPI.sh ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ../ljmd.mpi.x)
endif()

# does it create the file write.ppm?
#file(WRITE ${CMAKE_BINARY_DIR}/test_write.sh "rm -f *.ppm && ./canvas-draw && test -f white.ppm")
#add_test(NAME Writes COMMAND bash test_write.sh)
