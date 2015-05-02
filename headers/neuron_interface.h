// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#ifndef NEURODYNAMICS_NEURON_INTERFACE_H_
#define NEURODYNAMICS_NEURON_INTERFACE_H_

#pragma once

#include "connections_interface.h"

class NeuronInterface {
 public:
  NeuronInterface() {}
  virtual int process(int turns) = 0;
  virtual ~NeuronInterface() {}
};
#endif  // NEURODYNAMICS_NEURON_INTERFACE_H_
