cmake_minimum_required(VERSION 3.10)

#enable_testing()

#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

#enable_language(CXX)
#include(GTest)

#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression.sh
#	DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/regression.sh
#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression.sh
#	DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/regression.sh

#add_custom_command(
#    TARGET test POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/tests/regression.sh ${CMAKE_CURRENT_BINARY_DIR}/tests/regression.sh)

#test: utilities
add_executable(test_utilities.x test_utilities.cpp)
target_include_directories(test_utilities.x PRIVATE ../include)
target_link_libraries(test_utilities.x PRIVATE GTest::GTestMain GTest::GTest *ljmd_rt)

add_test(NAME Base COMMAND test_utilities.x)

#does the application run?
add_test(NAME Runs COMMAND echo "" > ../ljmd.x)

add_test(NAME Regression COMMAND regression.sh)

# does it create the file write.ppm?
#file(WRITE ${CMAKE_BINARY_DIR}/test_write.sh "rm -f *.ppm && ./canvas-draw && test -f white.ppm")
#add_test(NAME Writes COMMAND bash test_write.sh)
