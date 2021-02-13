cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.1 CACHE STRING "Floating point number for the current version")

project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    #set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
#message(STATUS "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

find_package(MPI REQUIRED)

#src files
list(APPEND src_files src/engine.c src/defs.c src/input_output.c src/utilities.c)
list(APPEND ccflags "-Wall" "-Wextra" "-DLJMD_VERSION=${LJMD_VERSION}")

add_library(ljmd_rt SHARED ${src_files})
target_include_directories(ljmd_rt PRIVATE include)
target_link_libraries(ljmd_rt PRIVATE m MPI::MPI_C)

add_executable(ljmd.x main.c)
target_include_directories(ljmd.x PRIVATE include)
target_compile_options(ljmd.x PRIVATE ${ccflags})
target_link_libraries(ljmd.x PRIVATE m ljmd_rt MPI::MPI_C)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()

    set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
    enable_language(CXX)
    include(GTest)

    add_subdirectory(tests)
endif()
