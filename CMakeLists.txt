cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.0 CACHE STRING "Floating point number for the current version")
project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

add_executable(ljmd.x src/ljmd.c)

target_compile_options(ljmd.x PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer)
target_compile_definitions(ljmd.x PRIVATE -DLJMD_VERSION=${LJMD_VERSION})
target_link_libraries(ljmd.x PRIVATE m)

set(ENABLE_TESTING on)

if(ENABLE_TESTING)
	enable_testing()

	set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
	
	enable_language(CXX)
	include(GTest)
	
	#add_custom_command(
    #    TARGET test POST_BUILD
    #    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/tests/regression.sh ${CMAKE_CURRENT_BINARY_DIR}/tests/regression.sh)
	
	#test: utilities
	add_executable(test_utilities.x tests/test_utilities.cpp)
	target_include_directories(test_utilities.x PRIVATE include)
	target_link_libraries(test_utilities.x PRIVATE GTest::GTestMain GTest::GTest)
	add_test(NAME Base COMMAND tests/test_utilities.x)

	#does the application run?
	add_test(NAME Runs COMMAND echo "" > ljmd.x)

	add_test(NAME Regression COMMAND tests/regression.sh)

	# does it create the file write.ppm?
	#file(WRITE ${CMAKE_BINARY_DIR}/test_write.sh "rm -f *.ppm && ./canvas-draw && test -f white.ppm")
	#add_test(NAME Writes COMMAND bash test_write.sh)
	
endif()
