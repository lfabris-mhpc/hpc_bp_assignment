cmake_minimum_required(VERSION 3.10)

set(LJMD_VERSION 1.0 CACHE STRING "Floating point number for the current version")
project(ljmd VERSION ${LJMD_VERSION} LANGUAGES C)

add_executable(ljmd-serial.x src/ljmd.c)

target_compile_options(ljmd-serial.x PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer)
target_compile_definitions(ljmd-serial.x PRIVATE -DLJMD_VERSION=${LJMD_VERSION})
target_link_libraries(ljmd-serial.x PRIVATE m)
