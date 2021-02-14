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
  
OpenMP branch (before merging with the Optimization branch): measured in a Linux VM running on a MacOS machine with Intel(R) Core(TM) i7-4510U CPU @ 2.00GHz, L3 cache: 4 MiB.  
- 1 OpenMP thread
argon_108:  speedup: 4.122; reference time:  39.199s  
argon_2916: speedup: 2.395; reference time: 966.978s  
- 2 OpenMP threads
argon_108:  speedup: 5.557; reference time:  39.199s  
argon_2916: speedup: 2.484; reference time: 966.978s
- 4 OpenMP threads
argon_108:  speedup: 5.526; reference time:  39.199s  
argon_2916: speedup: 2.501; reference time: 966.978s  
  
MPI branch (before merging with the Optimizations branch): measured on a machine with Intel(R) Core(TM) i7-7700HQ CPU @ 2.80GHz, L3 cache: 6 MiB.  
Different implementations of the force calculation were benchmarked with 4 MPI processes, with a Release build; reference time are from the base version  
- parallelization as seen in class  
argon_108:  speedup: 6.751; reference time:   7.904s  
argon_2916: speedup: 6.487; reference time: 152.732s  
- as before, but using non-blocking communications  
argon_108:  speedup: 6.422  
argon_2916: speedup: 6.137  
- as before, but each process is assigned different ranges to achieve balanced estimated workloads  
argon_108:  speedup: 6.117  
argon_2916: speedup: 5.809  
- naive reimplementation without symmetry  
argon_108:  speedup: 3.276  
argon_2916: speedup: 2.907  
- naive reimplementation with non-blocking communications, no symmetry  
argon_108:  speedup: 3.270  
argon_2916: speedup: 3.039  
- ring pattern overlapping communication and computation, no symmetry  
argon_108:  speedup: 3.255  
argon_2916: speedup: 2.914  
- ring pattern overlapping communication and computation, exploiting symmetry  
argon_108:  speedup: 3.748  
argon_2916: speedup: 3.805  
  
Main branch benchmarks after integration of the other branches: measured on a machine with Intel(R) Core(TM) i7-7700HQ CPU @ 2.80GHz, L3 cache: 6 MiB.  
Parallelization using a total of 4 processes/threads, with a Release build; reference time are from the base version  
- 4 MPI processes  
argon_108:  speedup: 17.769; reference time:   7.904s  
argon_2916: speedup:  8.211; reference time: 152.732s  
- 2 MPI processes, 2 OpenMP threads  
argon_108:  speedup: 17.793  
argon_2916: speedup:  8.259  
- 4 OpenMP threads  
argon_108:  speedup: 18.245  
argon_2916: speedup:  8.550  
  
The contributors are:  
- Mattia Mencagli as mattiamencagli  
- Orlenys Troconis as otroconi, student 
- Lorenzo Fabris as lfabris-mhpc  
