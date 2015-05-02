// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#ifndef NEURODYNAMICS_NEURONETWORK_H_
#define NEURODYNAMICS_NEURONETWORK_H_

#pragma once

#include "connections_interface.h"
#include "connections_sparse.h"
#include "connections_dense.h"

#include "neuron_interface.h"
#include "neuron_hodgin.h"
#ifdef NEURODYNAMICS_WITH_MPI
#include "connections_dense_mpi.h"
#include "neuron_hodgin_mpi.h"
#endif

class Neuronetwork {
 private:
  NeuronInterface* Neurons;
  ConnectionsInterface* Connections;
  bool initialized_outside = false;

 public:
  Neuronetwork(NeuronInterface* _Neurons, ConnectionsInterface* _Connections);

  Neuronetwork(int neuron_num);
  Neuronetwork(int neuron_num, int receptor_type);
  Neuronetwork(int neuron_num, int receptor_type,
                             int external_function);
  int Process(int turns);
  ~Neuronetwork();
};
#endif  // NEURODYNAMICS_NEURONETWORK_H_
