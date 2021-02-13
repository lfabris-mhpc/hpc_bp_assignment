This package contains simplified MD code with multi-threading and multi-process 
parallelization for simulating atoms with a Lennard-Jones potential.  

The project is built with CMake and provides test through the GoogleTest library.  

The examples directory contains 2 sets of example input decks
and the reference directory the corresponding outputs.  

Type:  
`$ mkdir build`  
`$ cd build`  
`$ cmake ..`  
optionally, select the build type with:  
`$ cmake -DCMAKE_BUILD_TYPE=<type> ..`  
with <type> one of Debug, Release, MinSizeRel, RelWithDebInfo  
to compile:  
`$ cmake --build .`  
to remove all compiled objects:  
`$ cmake --build . --target clean`  

Benchmarking (the speedup is measured by dividing the time obtained in the v.0.base by the time of the benchmarked branch):  
Optimization branch: measured on a machine with Intel(R) Core(TM) i5-8250U CPU @ 1.60GHz, L3 cache: 6 MiB.  
argon_108:  speedup : 3.578; reference time:  11.243s   
argon_2916: speedup : 2.560; reference time: 288.341s
  
OpenMP branch (before merging with the Optimization branch): measured on a machine with Intel(R) Core(TM) i7-4510U CPU @ 2.00GHz, L3 cache: 4 MiB.  
nthreads=1
argon_108:  speedup: 4.122; reference time:  39.199s  
argon_2916: speedup: 2.395; reference time: 966.978s  
nthreads=2
argon_108:  speedup: 5.557; reference time:  39.199s  
argon_2916: speedup: 2.484; reference time: 966.978s
nthreads=4
argon_108:  speedup: 5.526; reference time:  39.199s  
argon_2916: speedup: 2.501; reference time: 966.978s
The contributors are:  

- Mattia Mencagli as mattiamencagli  
- Orlenys Troconis as otroconi  
- Lorenzo Fabris as lfabris-mhpc  
