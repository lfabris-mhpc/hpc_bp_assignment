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

The contributors are:
- Mattia Mencagli as mattiamencagli
- Orlenys Troconis as otroconi
- Lorenzo Fabris as lfabris-mhpc
