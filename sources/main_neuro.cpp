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

#ifdef NEURODYNAMICS_WITH_MPI
  MPI_Init(&argc, &argv);
  int num_ranks;
  int my_rank;
  double start_t, init_finish_t, process_finish_t;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#ifdef NEURODYNAMICS_DEBUG
  std::cout << "Hello from rank " << my_rank << " out of " << num_ranks << std::endl;
#endif
  start_t = MPI_Wtime();
  ConnectionsDenseMPI *connections = new ConnectionsDenseMPI(my_rank, num_ranks, neuron_num, connection_probability,
                                                             AMPA_RECEPTOR);
  NeuronHodgkinMPI *neurons = new NeuronHodgkinMPI(neuron_num, connections, my_rank, num_ranks, AMPA_RECEPTOR,
                                                   I_EXTERNAL_RANDOM);
  Neuronetwork hh(neurons, connections);

  if (my_rank == 0) {
      std::cout << "Iterations             " << turns_num << std::endl;
      std::cout << "Amount of neurons      " << neuron_num << std::endl;
      std::cout << "Connection probability " << connection_probability << std::endl;
  }
  init_finish_t = MPI_Wtime();
#else
  clock_t start_t, init_finish_t, process_finish_t;
  start_t = clock();
  ConnectionsDense* connections = new ConnectionsDense(neuron_num, connection_probability, AMPA_RECEPTOR);
  NeuronHodgkin* neurons = new NeuronHodgkin(neuron_num, connections, AMPA_RECEPTOR, I_EXTERNAL_RANDOM);
  Neuronetwork hh(neurons, connections);
  init_finish_t = clock();
#endif

  hh.Process(turns_num); // main cycle

#ifdef NEURODYNAMICS_WITH_MPI
  process_finish_t = MPI_Wtime();
  double local_init_time = init_finish_t - start_t;
  double local_process_time = process_finish_t - init_finish_t;
  if (my_rank < 64) {
    std::cout << "Rank " << my_rank << ": Initilization took " << local_init_time << std::endl;
    std::cout << "Rank " << my_rank << ": Process took " << local_process_time << std::endl;
  }
  double lowest_time, highest_time, sum_of_time;
  MPI_Reduce(&local_init_time, &lowest_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_init_time, &sum_of_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_init_time, &highest_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (my_rank == 0) {
    std::cout << "Lowest init time " << lowest_time << std::endl;
    std::cout << "Average init time " << sum_of_time/num_ranks << std::endl;
    std::cout << "Highest init time " << highest_time << std::endl;
  }
  MPI_Reduce(&local_process_time, &lowest_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_process_time, &sum_of_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_process_time, &highest_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (my_rank == 0) {
    std::cout << "Lowest process time " << lowest_time << std::endl;
    std::cout << "Average process time " << sum_of_time/num_ranks << std::endl;
    std::cout << "Highest process time " << highest_time << std::endl;
  }
  MPI_Finalize();
#else
  process_finish_t = clock();
  std::cout << " Initilization took " << static_cast<double>(init_finish_t - start_t) / CLOCKS_PER_SEC << std::endl;
  std::cout << " Process took " << static_cast<double>(process_finish_t - init_finish_t) / CLOCKS_PER_SEC << std::endl;
#endif
  delete neurons;
  delete connections;
  return 0;
}
