cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.0 CACHE STRING "Floating point number for the current version")

project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

#src files
list(APPEND src_files src/engine.c src/defs.c src/input_output.c src/utilities.c)

#basic compiler flags
list(APPEND ccflags "-Wall" "-Wextra" "-DLJMD_VERSION=${LJMD_VERSION}")

#actual targets
message(STATUS "Setting basic serial targets")
add_library(ljmd_rt SHARED ${src_files})
target_include_directories(ljmd_rt PRIVATE include)
target_link_libraries(ljmd_rt PRIVATE m)

add_executable(ljmd.x main.c)
target_include_directories(ljmd.x PRIVATE include)
target_compile_options(ljmd.x PRIVATE ${ccflags})

target_link_libraries(ljmd.x PRIVATE m ljmd_rt)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
	
	set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
	enable_language(CXX)
	include(GTest)
	
    add_subdirectory(tests)
endif()
