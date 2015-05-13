// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#include "defines.h"
#include "neuronetwork.h"

Neuronetwork::Neuronetwork(NeuronInterface* _Neurons, ConnectionsInterface* _Connections) {
  Neurons = _Neurons;
  Connections = _Connections;
  initialized_outside = true;
}

Neuronetwork::Neuronetwork(int neuron_num) {
  Connections = new ConnectionsDense(neuron_num, CONNECTION_PROBABILITY);
  Neurons = new NeuronHodgkin(neuron_num, Connections, DEFAULT_RECEPTOR);
}

Neuronetwork::Neuronetwork(int neuron_num, int receptor_type) {
  Connections = new ConnectionsDense(neuron_num, CONNECTION_PROBABILITY,
                                     receptor_type);
  Neurons = new NeuronHodgkin(neuron_num, Connections, receptor_type);

}

Neuronetwork::Neuronetwork(int neuron_num, int receptor_type,
                           int external_function) {
  Connections = new ConnectionsDense(neuron_num, CONNECTION_PROBABILITY,
                                     receptor_type);
  Neurons = new NeuronHodgkin(neuron_num, Connections, receptor_type,
                              external_function);
}

int Neuronetwork::Process(int turns) {
  Neurons->process(turns);
  return 0;
}

Neuronetwork::~Neuronetwork() {
  if (!initialized_outside) {
    delete Neurons;
    delete Connections;
  }
}
