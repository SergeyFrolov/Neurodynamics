// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#ifndef NEURODYNAMICS_CONNECTIONS_DENSE_H_
#define NEURODYNAMICS_CONNECTIONS_DENSE_H_

#pragma once

#include "connections_interface.h"

class ConnectionsDense : public ConnectionsInterface {
 protected:
  double* connections;
  /*
  0 1 1 1
  0 0 0 1
  1 0 0 1
  0 0 0 0
  1st  neuron sends to  everyone.
  last neuron gets from everyone.
  diagonal is always zero
  (or maybe not always, if it makes biological sense)
  */
  void GetNeuronIdsToSend(int neuron_id, std::vector<int>* send_ids);
  void GetNeuronsToRecv(int neuron_id, std::vector<SynapticConnection>* recv_connections);
 public:
  // _size is a size of single dimension, resulting size would be size*size
  // density is a double between 0 and 1
  ConnectionsDense(int _size, int graph_type = GRAPH_DEFAULT,
                   double density = GRAPH_RAND_P_CONNECTION,
                   int receptor_type = AMPA_RECEPTOR);
  ConnectionsDense() {}
  virtual void GetAllNeuronIdsToSend(std::vector<std::vector<int>>* send_ids);
  virtual void GetAllNeuronsToRecv(std::vector<std::vector<SynapticConnection>>* recv_connections);
  void print();
 virtual ~ConnectionsDense();
};
#endif  // NEURODYNAMICS_CONNECTIONS_DENSE_H_
