cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.0 CACHE STRING "Floating point number for the current version")

project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    #set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
#message(STATUS "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

list(APPEND src_files src/engine.c src/defs.c src/input_output.c src/utilities.c)

list(APPEND ccflags "-Wall" "-Wextra" "-DLJMD_VERSION=${LJMD_VERSION}")
list(APPEND ccflags_copt "-g" "-O3" "-march=native" "-DNDEBUG" "-ffast-math" "-fomit-frame-pointer" "-fexpensive-optimizations")

#force selection
if(FORCE)
	list(APPEND ccflags "-DFORCE=${FORCE}")
endif()

add_library(ljmd_rt SHARED ${src_files})

target_include_directories(ljmd_rt PRIVATE include)
target_compile_options(ljmd_rt PRIVATE ${ccflags} ${ccflags_copt})
target_link_libraries(ljmd_rt PRIVATE MPI::MPI_C OpenMP::OpenMP_C m)

add_executable(ljmd.x main.c)
target_include_directories(ljmd.x PRIVATE include)

target_compile_options(ljmd.x PRIVATE ${ccflags} ${ccflags_copt})
target_link_libraries(ljmd.x PRIVATE m ljmd_rt MPI::MPI_C OpenMP::OpenMP_C)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()

    set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
    enable_language(CXX)
    include(GTest)

    add_subdirectory(tests)
endif()
