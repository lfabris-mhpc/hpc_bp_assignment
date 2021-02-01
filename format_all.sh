#!/bin/bash

clang-format -i --style=file main.c
clang-format -i --style=file src/*.c
clang-format -i --style=file include/*.h

clang-format -i --style=file tests/*.cpp
#for the future
#clang-format -i --style=file tests/*.hpp
#clang-format -i --style=file tests/*.h
