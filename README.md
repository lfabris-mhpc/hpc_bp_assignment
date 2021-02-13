This package contains simplified MD code with multi-threading and multi-process 
parallelization for simulating atoms with a Lennard-Jones potential.

The project is built with CMake and provides test through the GoogleTest library.

The examples directory contains 2 sets of example input decks
and the reference directory the corresponding outputs.

Type: mkdir build
$ cd build
$ cmake ..
optionally, select the build type with
$ cmake .. -DCMAKE_BUILD_TYPE=<type>
with <type> one of Debug, Release, MinSizeRel, RelWithDebInfo
to compile
$ cmake --build .
to remove all compiled objects
$ cmake --build . --target clean
  
Benchmarking (the speedup is measured by dividing the time obtained in the v.0.base by the time of the benchmarked branch):
- Optimization branch, argon_108:  Speedup = 3.578
- Optimization branch, argon_2916: Speedup = 2.560

The contributors are:
- Mattia Mencagli as mattiamencagli
- Orlenys Troconis as otroconi
- Lorenzo Fabris as lfabris-mhpc
