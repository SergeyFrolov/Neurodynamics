// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "neuron_hodgin.h"

NeuronHodgkin::NeuronHodgkin(unsigned int first_neuron, unsigned int last_neuron,
                             ConnectionsInterface* Connections,
                             unsigned int _receptor_type, unsigned int _external_function) {
  this->receptor_type = _receptor_type;
  this->neuron_num = last_neuron - first_neuron + 1;

  V_old = new double[neuron_num];
  V_new = new double[neuron_num];
  m     = new double[neuron_num];
  h     = new double[neuron_num];
  n     = new double[neuron_num];
  I_ext = new double[neuron_num];
  I_syn = new double[neuron_num];
  send_Sact = new double[neuron_num];

  recv_connections.resize(neuron_num);
  send_ids.resize(neuron_num);
  Connections->GetAllNeuronIdsToSend(&(send_ids));
  Connections->GetAllNeuronsToRecv(&(recv_connections));

  presynaptic_neurons = new remote_neuron*[neuron_num];
  for (unsigned int i = first_neuron; i < first_neuron + neuron_num; i++) {
    presynaptic_neurons[i] = new remote_neuron[recv_connections[i].size()];
  }

  switch (_external_function) {
    case I_EXTERNAL_NULL:
      for (unsigned int i = first_neuron; i < first_neuron + neuron_num; i++) {
        I_ext[i] = 0.0;
      }
      break;

    case I_EXTERNAL_RANDOM:
      for (unsigned int i = first_neuron; i < first_neuron + neuron_num; i++) {
        I_ext[i] = (static_cast<double>(rand()) / RAND_MAX)
                   * I_EXTERNAL_RANDOM_MAX_VALUE;
      }
      break;

    default:
      throw;
  }

  for (unsigned int i = first_neuron; i < first_neuron + neuron_num; i++) {
    V_old[i] = V_new[i] = V_REST;
    I_syn[i] = 0;

    m[i] = 0.0526;
    h[i] = 0.598;
    n[i] = 0.317;
    send_Sact[i] = 0.0001;
  }
}

NeuronHodgkin::NeuronHodgkin(unsigned int _neuron_num,
                             ConnectionsInterface* Connections,
                             unsigned int _receptor_type, unsigned int _external_function) :
  NeuronHodgkin(0, _neuron_num - 1, Connections, _receptor_type, _external_function) {
}


int NeuronHodgkin::process(int turns) {
  double* tmp_ptr_V;
  print(-1, OUTPUT_INIT, 3);

  for (int cur_turn = 0; cur_turn < turns; cur_turn++) {
    for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
      send_Sact[current_neuron] = CalcPostSActivation(send_Sact[current_neuron]);
      SendRecvPresynaptic();
    }
    for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
      m[current_neuron]     = RungeKutta4(&NeuronHodgkin::_CalcDerivativeFm,
                                          m[current_neuron], DELTAT);
      h[current_neuron]     = RungeKutta4(&NeuronHodgkin::_CalcDerivativeFh,
                                          h[current_neuron], DELTAT);
      n[current_neuron]     = RungeKutta4(&NeuronHodgkin::_CalcDerivativeFn,
                                          n[current_neuron], DELTAT);
    }
    for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
      I_syn[current_neuron] = CalcSynapticCurrent();
      V_new[current_neuron] = RungeKutta4(&NeuronHodgkin::_CalcDerivativeFv,
                                          V_old[current_neuron], DELTAT);
    }

    tmp_ptr_V = V_new;
    V_new = V_old;
    V_old = tmp_ptr_V;
    print(cur_turn);
  }

  print(-1, OUTPUT_FINAL, 3);
  return 0;
}

void NeuronHodgkin::SendRecvPresynaptic() {
  for (unsigned int iter_neuron = 0; iter_neuron < recv_connections[current_neuron].size(); iter_neuron++) {
    presynaptic_neurons[current_neuron][iter_neuron].voltage =
      V_old[recv_connections[current_neuron][iter_neuron].id];
    presynaptic_neurons[current_neuron][iter_neuron].Sact_generated =
      send_Sact[recv_connections[current_neuron][iter_neuron].id];
  }
}
/*
Returns synaptic current, generated on all connections.
  recv_connections[current_neuron][iter_neuron].g  - max conductance
  --/--.Sact_generated  - activation variable. Percentage of neurotransmitter docked on the
                          postsynaptic cell relative to the maximum that can dock.
  RECEPTOR_TYPE is from defines.h
*/
double NeuronHodgkin::CalcSynapticCurrent() {
  double Esyn;

  switch (receptor_type) {
    case AMPA_RECEPTOR:
      Esyn = AMPA_E;
      break;

    case NMDA_RECEPTOR:
    case GABA_A_RECEPTOR:
    case GABA_B_RECEPTOR:
    default:
      throw;
  }

  double result = 0;
  for (unsigned int iter_neuron = 0; iter_neuron < recv_connections[current_neuron].size(); iter_neuron++) {
    result += recv_connections[current_neuron][iter_neuron].g *
              presynaptic_neurons[current_neuron][iter_neuron].Sact_generated *
              (V_old[current_neuron] - Esyn);
  }
  return result;
}

double NeuronHodgkin::_CalcDerivativeFv(double V) {
  return -(_I_Na(V) + _I_K(V) + _I_L(V) + I_syn[current_neuron] -
           I_ext[current_neuron]) / C_M;
}

/* REGION: Currents */
// sodium(Na)
inline double NeuronHodgkin::_I_Na(double V) {
  return G_NA_MAX
         * m[current_neuron] * m[current_neuron] * m[current_neuron]
         * h[current_neuron]
         * (V - V_NA_REST);
}

inline double NeuronHodgkin::_I_K(double V) {
  return G_K_MAX
         * n[current_neuron] * n[current_neuron] * n[current_neuron] * n[current_neuron]
         * (V - V_K_REST);
}

inline double NeuronHodgkin::_I_L(double V) {
  return G_L
         * (V - V_L_REST);
}
/* END OF REGION: Currents */

/* REGION: Gate variables */
#define devV (V_old[current_neuron] - V_REST)
// sodium(Na) channel goes through 2 gates: m and h.
double NeuronHodgkin::_CalcDerivativeFm(double m_var) {
  double alpha_m = (2.5 - 0.1 * devV) / (std::exp(2.5 - 0.1 * devV) - 1);
  double beta_m = 4 / std::exp(devV / 18);
  return alpha_m * (1 - m_var) - beta_m * m_var;
}

double NeuronHodgkin::_CalcDerivativeFh(double h_var) {
  double alpha_h = 0.07 / std::exp(0.05 * devV);
  double beta_h = 1 / (std::exp(3 - 0.1 * devV) + 1);
  return alpha_h * (1 - h_var) - beta_h * h_var;
}

// potassium(K) channel goes through 1 gate: n.
double NeuronHodgkin::_CalcDerivativeFn(double n_var) {
  double alpha_n = (0.1 - 0.01 * devV) / (std::exp(1 - 0.1 * devV) - 1);
  double beta_n = 0.125 / std::exp(0.0125 * devV);
  return alpha_n * (1 - n_var) - beta_n * n_var;
}
/* END OF REGION: Gate variables */

/* REGION: S activation */
double NeuronHodgkin::_CalcAmpaPotential(double S) {
  double Tconc = RECEPTOR_Tmax /
                 (1. + std::exp(-(V_old[current_neuron] - RECEPTOR_Vp) / RECEPTOR_Kp));
  return AMPA_ALPHA * Tconc * (1. - S) - AMPA_BETA * S;
}
double NeuronHodgkin::_CalcNdmaPotential(double S) {
  throw;
  return S;
}
double NeuronHodgkin::_CalcGabaAPotential(double S) {
  throw;
  return S;
}
double NeuronHodgkin::_CalcGabaBPotential(double S) {
  throw;
  return S;
}

/*
Activation variable calculation. Returns next Sact for postsynaptic neuron.
Sact  - previous Sact. Percentage of neurotransmitter docked on the
postsynaptic cell relative to the maximum that can dock.

Warning: the fact that Sact depends only on presynaptic voltage(and previous Sact)
gives opportunity for tricky optimization: we could calculate Sact of all postsynaptic
neurons on presynaptic neuron once and send it.
*/
double NeuronHodgkin::CalcPostSActivation(double Sact) {
  double (NeuronHodgkin::*PCalcPotential)(double);

  switch (receptor_type) {
    case AMPA_RECEPTOR:
      PCalcPotential = &NeuronHodgkin::_CalcAmpaPotential;
      break;

    case NMDA_RECEPTOR:
      PCalcPotential = &NeuronHodgkin::_CalcNdmaPotential;
      break;

    case GABA_A_RECEPTOR:
      PCalcPotential = &NeuronHodgkin::_CalcGabaAPotential;
      break;

    case GABA_B_RECEPTOR:
      PCalcPotential = &NeuronHodgkin::_CalcGabaBPotential;
      break;

    default:
      throw;
  }
  return RungeKutta4(PCalcPotential, Sact, DELTAT);
}
/* END OF REGION: S activation */

double NeuronHodgkin::RungeKutta4(double (NeuronHodgkin::*f)(double),
                                  double x, double dx) {
  double
  k1 = dx * (this->*f)(x),
  k2 = dx * (this->*f)(x + k1 / 2),
  k3 = dx * (this->*f)(x + k2 / 2),
  k4 = dx * (this->*f)(x + k3);
  return x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}


void NeuronHodgkin::print(int step, std::string name, int process_level) {
  if (process_level < 1)
  { return; }
  std::streambuf *coutbuf;
  std::ofstream fmid(name, std::ios_base::app);
  coutbuf = std::cout.rdbuf(); // save old buf
  std::cout.rdbuf(fmid.rdbuf()); // redirect std::cout

  if (step <= 0) {
    std::ofstream fmid(name);
    fmid.close();
  }
  if (process_level == 1) {

    for (unsigned int i = 0; i < neuron_num; i++)
      if ((V_old[i] > V_new[i]) && (V_old[i] > OUTPUT_PROCESS_PEAK_THRESHOLD)) {

        std::cout << "Step " << step << ": neuron " << i <<
                  " has " << V_old[i] << " mV" << std::endl;
      }

  }

  if (process_level == 2) {
    for (unsigned int i = 0; i < neuron_num; i++) {
      std::cout << V_new[i] << " ";
    }

    std::cout << std::endl;
  }
  if (process_level >= 3) {
    if (step >= 0)
    { std::cout << "Step " << step << std::endl; }
    std::cout << "V: ";

    for (unsigned int i = 0; i < neuron_num; i++) {
      std::cout << V_new[i] << " ";
    }

    std::cout << std::endl;
    std::cout << "I_ext: ";

    for (unsigned int i = 0; i < neuron_num; i++) {
      std::cout << I_ext[i] << " ";
    }

    std::cout << std::endl;
    std::cout << "I_syn: ";

    for (unsigned int i = 0; i < neuron_num; i++) {
      std::cout << I_syn[i] << " ";
    }

    std::cout << std::endl;
    std::cout << "m: ";

    for (unsigned int i = 0; i < neuron_num; i++) {
      std::cout << m[i] << " ";
    }

    std::cout << std::endl;
    std::cout << "h: ";

    for (unsigned int i = 0; i < neuron_num; i++) {
      std::cout << h[i] << " ";
    }

    std::cout << std::endl;
    std::cout << "n: ";

    for (unsigned int i = 0; i < neuron_num; i++) {
      std::cout << n[i] << " ";
    }

    std::cout << std::endl;
  }
  if (process_level >= 4) {
    std::cout << "send_S_act:\n";
    for (unsigned int i = 0; i < neuron_num; i++) {
      std::cout << send_Sact[i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout.rdbuf(coutbuf); // reset to standard output again
}

NeuronHodgkin::~NeuronHodgkin() {
  delete[] V_old;
  delete[] V_new;
  delete[] m;
  delete[] h;
  delete[] n;
  delete[] I_ext;
  delete[] I_syn;
  for (unsigned int i = 0; i < neuron_num; i++) {
    if (presynaptic_neurons[i] != NULL)
      delete[] presynaptic_neurons[i];
  }
  delete[] presynaptic_neurons;
  delete[] send_Sact;
}
