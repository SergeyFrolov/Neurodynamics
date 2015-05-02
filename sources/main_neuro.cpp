// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#include <stdlib.h>
#include <time.h>
#include <map>
#include <list>
#include <unistd.h>
#include <iostream>
#include <limits>

#include "defines.h"
#include "neuronetwork.h"

#ifdef NEURODYNAMICS_WITH_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[]) {
  dup2(STDOUT_FILENO, STDERR_FILENO);// redirecting stderr to stdout
#if RANDOMIZE_SEED
  srand(unsigned(time(NULL)));
#else
  srand(0);
#endif

  int turns_num = (argc >= 2) ? atoi(argv[1]) : PROCESS_TURNS;
  int neuron_num = (argc >= 3) ? atoi(argv[2]) : NEURON_NUM;
  double connection_probability = (argc >= 4) ? atof(argv[3]) : CONNECTION_PROBABILITY;

  clock_t start_t, init_finish_t, process_finish_t;
  start_t = clock();

#ifdef NEURODYNAMICS_WITH_MPI
  MPI_Init(&argc, &argv);
  int num_ranks;
  int my_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#ifdef NEURODYNAMICS_DEBUG
  std::cout << "Hello from rank " << my_rank << " out of " << num_ranks << std::endl;
#endif
  ConnectionsDenseMPI *connections = new ConnectionsDenseMPI(my_rank, num_ranks, neuron_num, connection_probability,
                                                             AMPA_RECEPTOR);
  NeuronHodgkinMPI *neurons = new NeuronHodgkinMPI(neuron_num, connections, my_rank, num_ranks, AMPA_RECEPTOR,
                                                   I_EXTERNAL_RANDOM);
  Neuronetwork hh(neurons, connections);
  if (my_rank == 0) {
      std::cout << "Iterations" << turns_num << std::endl;
      std::cout << "Amount of neurons " << neuron_num << std::endl;
      std::cout << "Connection probability " << connection_probability << std::endl;
  }
#else
  ConnectionsDense* connections = new ConnectionsDense(neuron_num, connection_probability, AMPA_RECEPTOR);
  NeuronHodgkin* neurons = new NeuronHodgkin(neuron_num, connections, AMPA_RECEPTOR, I_EXTERNAL_RANDOM);
  Neuronetwork hh(neurons, connections);
#endif

  init_finish_t = clock();
  hh.Process(turns_num);
  process_finish_t = clock();

  delete neurons;
  std::cout << "Rank " << my_rank << ": Deleted neurons." << std::endl;
  delete connections;
  std::cout << "Rank " << my_rank << ": Deleted connections." << std::endl;
#ifdef NEURODYNAMICS_WITH_MPI
  std::cout << "Rank " << my_rank << ": Initilization took " <<
  static_cast<double>(init_finish_t - start_t) / CLOCKS_PER_SEC << std::endl;
  std::cout << "Rank " << my_rank << ": Process took " <<
  static_cast<double>(process_finish_t - init_finish_t) / CLOCKS_PER_SEC << std::endl;
  MPI_Finalize();
#else
  std::cout << " Initilization took " << static_cast<double>(init_finish_t - start_t) / CLOCKS_PER_SEC << std::endl;
  std::cout << " Process took " << static_cast<double>(process_finish_t - init_finish_t) / CLOCKS_PER_SEC << std::endl;
#endif
  return 0;
}
