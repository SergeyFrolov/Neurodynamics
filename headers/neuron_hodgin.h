// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#ifndef NEURODYNAMICS_NEURON_HODGIN_H_
#define NEURODYNAMICS_NEURON_HODGIN_H_

#pragma once

#include <vector>
#include <string>

#include "defines.h"
#include "neuron_interface.h"
#include "connections_interface.h"

struct remote_neuron {
  double voltage;
  double Sact_generated; // actually postsynaptic Sact. Calculted on presyn for optimization.
};

class NeuronHodgkin : public NeuronInterface {
 protected:
  unsigned int neuron_num;
  double *V_old, *V_new;
  double *m, *h, *n;
  double *I_ext;
  double *I_syn;
  unsigned int receptor_type;

  std::vector<std::vector<SynapticConnection>> recv_connections;
  remote_neuron** presynaptic_neurons;

  std::vector<std::vector<int>> send_ids;
  double* send_Sact; // Sact of neurons i'th neuron sends to

  unsigned int current_neuron;
  // since we want Fv, Fm, Fh, Fn to be functions of one variable
  // to be solvable with Runge-Kutta, we have to have a state.

 public:
  NeuronHodgkin() {}
  NeuronHodgkin(unsigned int first_neuron, unsigned int last_neuron,
                ConnectionsInterface* Connections,
                unsigned int _receptor_type,
                unsigned int _external_function = I_EXTERNAL_NULL);
  NeuronHodgkin(unsigned int _neuron_num, ConnectionsInterface* Connections,
                unsigned int _receptor_type,
                unsigned int _external_function = I_EXTERNAL_NULL);

  virtual int process(int turns);
  virtual ~NeuronHodgkin();

 protected:
  double _CalcDerivativeFv(double V);

  double _I_Na(double V);
  double _I_K(double V);
  double _I_L(double V);

  double _CalcDerivativeFm(double m);
  double _CalcDerivativeFh(double h);
  double _CalcDerivativeFn(double n);

  virtual void SendRecvPresynaptic();
  double CalcSynapticCurrent();
  double CalcPostSActivation(double Sact);
  double _CalcAmpaPotential(double S);
  double _CalcNdmaPotential(double S);
  double _CalcGabaAPotential(double S);
  double _CalcGabaBPotential(double S);

  double RungeKutta4(double(NeuronHodgkin::*f)(double), double x, double dx);

  virtual void print(int step, std::string name = OUTPUT_PROCESS, int process_level = OUTPUT_PROCESS_LEVEL);
};
#endif  // NEURODYNAMICS_NEURON_HODGIN_H_
