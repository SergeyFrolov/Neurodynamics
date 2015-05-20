// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#ifndef NEURODYNAMICS_CONNECTIONS_DENSE_MPI_H_
#define NEURODYNAMICS_CONNECTIONS_DENSE_MPI_H_

#pragma once

#include "connections_dense.h"
#include <mpi.h>

class ConnectionsDenseMPI : public ConnectionsDense {
 private:
  int my_rank;
  int num_ranks;
  int n_per_rank;
  int n_per_your_rank;

  MPI_Datatype MPI_REMOTE_NEURON;
 public:
  // _size is a size of single dimension, resulting size would be size*size
  // density is a double between 0 and 1
  ConnectionsDenseMPI(int _my_rank, int _num_ranks, int _size,
                      int graph_type = GRAPH_DEFAULT, double density = GRAPH_RAND_P_CONNECTION,
                      int receptor_type = AMPA_RECEPTOR);
  void GetAllNeuronIdsToSend(std::vector<std::vector<int>>* send_ids) override;
  void GetAllNeuronsToRecv(std::vector<std::vector<SynapticConnection>>* recv_connections) override;
  ~ConnectionsDenseMPI() {}
};
#endif  // NEURODYNAMICS_CONNECTIONS_DENSE_MPI_H_
