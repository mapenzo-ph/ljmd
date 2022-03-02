This package contains simplified MD code with multi-threading
parallelization for simulating atoms with a Lennard-Jones potential.


The examples directory contains 3 sets of example input decks
and the reference directory the corresponding outputs.


To build using CMake type:
- mkdir build && cd build -> to create build folder
- cmake .. 		  -> to generate build files
- cmake --build .	  -> to build in the folder

To build with testing enabled :
- mkdir build && cd build
- cmake -DENABLE_TESTING .. for macos:  cmake -DENABLE_TESTING=1 -DCMAKE_CXX_FLAGS="-std=c++11" ..
- cmake --build .
- make test



The individual parts of this project is done by the following members of the team:

- Optimisation of force computation: Zainab: zainabnazari
- MPI parallelization: Giulio: giumal
- OpenMP parallelization: Alessandro: mapenzo-ph
