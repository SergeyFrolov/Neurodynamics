// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#include <assert.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

#include "connections_dense_mpi.h"

ConnectionsDenseMPI::ConnectionsDenseMPI(int _my_rank, int _num_ranks,
                                         int _size,  int graph_type, double _density,
                                         int receptor_type) {

  assert(_size > 0);
  n_per_rank = (_size - 1) / _num_ranks + 1;
  n_per_your_rank = ((_my_rank == _num_ranks - 1) && (_size % _num_ranks != 0)) ?
                    _size % n_per_rank :
                    n_per_rank;
  my_rank = _my_rank;
  num_ranks = _num_ranks;

  assert(_density <= 1);
  assert(_density >= 0);
  density = _density;
  size = _size;

  if (_my_rank == 0) {
    unsigned int edges = 0;
    double gMin, gMax;
    switch (receptor_type) {
      case AMPA_RECEPTOR:
        gMin = gAMPA_MIN;
        gMax = gAMPA_MAX;
        break;
      case NMDA_RECEPTOR:
        gMin = gNMDA_MIN;
        gMax = gNMDA_MAX;
        break;
      case GABA_A_RECEPTOR:
        gMin = gGABAA_MIN;
        gMax = gGABAA_MAX;
        break;
      case GABA_B_RECEPTOR:
        gMin = gGABAB_MIN;
        gMax = gGABAB_MAX;
        break;
      default:
        throw;
    }
    connections = new double[size * size];
    switch (graph_type) {
        case GRAPH_RAND:
            for (unsigned int i = 0; i < size; i++)
                for (unsigned int j = 0; j < size; j++)
                    if ((i == j) || ((static_cast<double>(rand()) / RAND_MAX) >= density)) {
                        connections[ind(j, i, size)] = 0.0;
                    } else {
                        connections[ind(j, i, size)] = gMin + (gMax - gMin) * static_cast<double>(rand()) / RAND_MAX;
                        edges++;
                    }
            break;
        case GRAPH_RING:
            for (unsigned int i = 0; i < size; i++)
                for (unsigned int j = 0; j < size; j++)
                    if ((i + 1 == j) || (( i == size - 1 ) && (j == 0))) {
                        connections[ind(j, i, size)] = gMin + (gMax - gMin) * static_cast<double>(rand()) / RAND_MAX;
                        edges++;
                    } else {
                        connections[ind(j, i, size)] = 0.0;
                    }
            break;
        default:
            throw;
    }
#ifdef OUTPUT_VERBOSE
  print();
#endif
    density = (size <= 1) ? 0 : edges / (size * (size - 1));
  }
  else {
    connections = NULL;
  }
  MPI_Bcast(&density, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_REMOTE_NEURON);
  MPI_Type_commit(&MPI_REMOTE_NEURON);
}

void ConnectionsDenseMPI::GetAllNeuronIdsToSend(std::vector<std::vector<int>>* send_ids) {
  assert(send_ids->size() == this->n_per_your_rank);
  if (my_rank == 0) {
    std::vector<int> tmp_vector;
    int tmp_vector_size;
    for (int iter_n = 0; iter_n < size; iter_n++) {
      if (iter_n < n_per_your_rank) {
        GetNeuronIdsToSend(iter_n, &((*send_ids)[iter_n]));
      } else {
        GetNeuronIdsToSend(iter_n, &tmp_vector);
        tmp_vector_size = tmp_vector.size();
        MPI_Send(&(tmp_vector_size), 1, MPI_INT, iter_n / n_per_rank, size + 1 + iter_n, MPI_COMM_WORLD);
        MPI_Send(tmp_vector.data(), tmp_vector_size, MPI_INT, iter_n / n_per_rank, iter_n, MPI_COMM_WORLD);
        tmp_vector.clear();
      }
    }
  } else {
    int vector_size;
    for (int iter_n = 0; iter_n < n_per_your_rank; iter_n++) {
      MPI_Recv(&(vector_size), 1, MPI_INT, 0, size + 1 + my_rank * n_per_rank + iter_n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      (*send_ids)[iter_n].resize(vector_size);
      MPI_Recv((*send_ids)[iter_n].data(), vector_size, MPI_INT, 0, my_rank * n_per_rank + iter_n, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }
}

void ConnectionsDenseMPI::GetAllNeuronsToRecv(std::vector<std::vector<SynapticConnection>>* recv_connections) {
  assert(recv_connections->size() == this->n_per_your_rank);
  if (my_rank == 0) {
    std::vector<SynapticConnection> tmp_vector;
    int tmp_vector_size;
    for (int iter_n = 0; iter_n < size; iter_n++) {
      if (iter_n < n_per_your_rank) {
        GetNeuronsToRecv(iter_n, &((*recv_connections)[iter_n]));
      } else {
        GetNeuronsToRecv(iter_n, &tmp_vector);
        tmp_vector_size = tmp_vector.size();
        MPI_Send(&(tmp_vector_size), 1, MPI_INT, iter_n / n_per_rank, size + 1 + iter_n, MPI_COMM_WORLD);
        MPI_Send(tmp_vector.data(), tmp_vector_size, MPI_REMOTE_NEURON, iter_n / n_per_rank, iter_n, MPI_COMM_WORLD);
        tmp_vector.clear();
      }
    }
  } else {
    int vector_size;
    for (int iter_n = 0; iter_n < n_per_your_rank; iter_n++) {
      MPI_Recv(&(vector_size), 1, MPI_INT, 0, size + 1 + my_rank * n_per_rank + iter_n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      (*recv_connections)[iter_n].resize(vector_size);
#ifdef NEURODYNAMICS_DEBUG
      std::cout << "Rank[" << my_rank << "] iter[" << iter_n <<
              "] receiving recv_connections of"  <<
              " size=" << vector_size <<
              " with tag " << my_rank * n_per_rank + iter_n <<
              ". Address of vector: "<< &((*recv_connections)[iter_n]) <<
              ". Address of data: "<< (*recv_connections)[iter_n].data() <<
              "." << std::endl;
#endif
      MPI_Recv((*recv_connections)[iter_n].data(), vector_size, MPI_REMOTE_NEURON, 0, my_rank * n_per_rank + iter_n,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
}
