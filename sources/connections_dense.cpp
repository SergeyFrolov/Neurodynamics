// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "defines.h"
#include "connections_dense.h"

ConnectionsDense::ConnectionsDense(int _size, double _density,
                                   int receptor_type) {

  assert(_size > 0);
  unsigned int edges = 0;

  assert(_density <= 1);
  assert(_density > 0);
  density = _density;

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
  size = _size;
  connections = new double[size * size];

  for (unsigned int i = 0; i < size; i++)
    for (unsigned int j = 0; j < size; j++)
      if ((i == j) || (rand() / RAND_MAX > density)) {
        connections[ind(j, i, size)] = 0.0;
      } else {
        connections[ind(j, i, size)] = gMin + (gMax - gMin) * static_cast<double>(rand()) / RAND_MAX;
        edges++;
      }

  print();

  density = (size == 1) ? 0 : edges / (size * (size - 1));
}

void ConnectionsDense::GetNeuronIdsToSend(int neuron_id, std::vector<int>* send_ids) {
  assert(send_ids->empty());
  for (unsigned int i = 0; i < size; i++) {
    if (connections[ind(i, neuron_id, size)] > CONNECT_EPS) {
      send_ids->push_back(i);
    }
  }
}

void ConnectionsDense::GetNeuronsToRecv(int neuron_id, std::vector<SynapticConnection>* recv_connections) {
  assert(recv_connections->empty());
  for (unsigned int i = 0; i < size; i++) {
    if (connections[ind(neuron_id, i, size)] > CONNECT_EPS) {
      recv_connections->push_back({ i, connections[ind(neuron_id, i, size)] });
    }
  }
}
void ConnectionsDense::GetAllNeuronIdsToSend(std::vector<std::vector<int>>* send_ids) {
  assert(send_ids->size() == size);
  for (unsigned int iter_n = 0; iter_n < size; iter_n++) {
    GetNeuronIdsToSend(iter_n, &(send_ids->operator[](iter_n)));
  }
}
void ConnectionsDense::GetAllNeuronsToRecv(std::vector<std::vector<SynapticConnection>>* recv_connections) {
  assert(recv_connections->size() == size);
  for (unsigned int iter_n = 0; iter_n < size; iter_n++) {
    GetNeuronsToRecv(iter_n, &(recv_connections->operator[](iter_n)));
  }
}

void ConnectionsDense::print() {
#ifdef OUTPUT_TOPOLOGY
  std::ofstream topoout(OUTPUT_TOPOLOGY);
  std::streambuf *coutbuf = std::cout.rdbuf();  // save old buf
  std::cout.rdbuf(topoout.rdbuf());  // redirect std::cout
#endif
  for (unsigned int i = 0; i < size; i++) {
    for (unsigned int j = 0; j < size; j++)
      std::cout << std::setw(6) << std::setprecision(3) <<
                connections[ind(j, i, size)] << " ";
    std::cout << std::endl;
  }
#ifdef OUTPUT_TOPOLOGY
  std::cout.rdbuf(coutbuf);  // reset to standard output again
#endif
}

ConnectionsDense::~ConnectionsDense() {
  if (connections != NULL)
    delete[] connections;
}
