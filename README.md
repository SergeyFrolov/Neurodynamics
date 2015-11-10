# Neurodynamics
Implementation of Hodgkin–Huxley model using Message Passing Interface.
## Build
CMakeFile.txt's for CMake are included. You may not use it, but, in any case, you will need c++ and/or MPI compilers(you could compile non-MPI version).   
## Usage   
This software could be used to simulate the dynamics of neurons(Hodgkin–Huxley model).   
Starting contitions are randomized, as I was more interested in impementing HPC part of it, but they may be easily initilaized in class, inherited from ConnectionsInterface or ConnectionsSparse(or you could just create another constructor. Essensially, any element > 0 in ConnectionsSparse is a directed connection between neurons with "power" equal to a value in connections[i][j].)   
NeuronHodgin is non-MPI version of class, handling dynamics.   
NeuronHodginMPI has 2 algorithms: p2p and collective. Almost always it is better to use collective.   
Scalabilty tested up to 200+ ranks.   
