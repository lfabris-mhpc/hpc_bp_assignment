cmake_minimum_required(VERSION 3.10)

#regression tests
file(COPY ../examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../reference DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY regression.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME Regression-test-serial COMMAND regression.sh 1)
add_test(NAME Regression-test-multip COMMAND regression.sh 4)

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
