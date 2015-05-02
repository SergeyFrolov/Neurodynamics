// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#ifndef NEURODYNAMICS_CONNECTIONS_INTERFACE_H_
#define NEURODYNAMICS_CONNECTIONS_INTERFACE_H_

#pragma once

#include <vector>
#include "defines.h"

struct SynapticConnection {
  unsigned int id;
  double g; // max conductance
};

class ConnectionsInterface {
 protected:
  unsigned int size;
  double density;
 public:
  ConnectionsInterface() {}

  virtual void GetNeuronIdsToSend(int neuron_id, std::vector<int>* send_ids) = 0;
  virtual void GetNeuronsToRecv(int neuron_id, std::vector<SynapticConnection>* recv_connections) = 0;
  virtual void GetAllNeuronIdsToSend(std::vector<std::vector<int>>* send_ids) = 0;
  virtual void GetAllNeuronsToRecv(std::vector<std::vector<SynapticConnection>>* recv_connections) = 0;
  virtual ~ConnectionsInterface() {}
  inline double GetDensity() { return density; }
};
#endif  // NEURODYNAMICS_CONNECTIONS_INTERFACE_H_
