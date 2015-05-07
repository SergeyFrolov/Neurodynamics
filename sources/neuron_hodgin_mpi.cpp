// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#include <stdlib.h>
#include <iostream>

#include "neuron_hodgin_mpi.h"
#include "mpi.h"

NeuronHodgkinMPI::NeuronHodgkinMPI(unsigned int _neuron_num, ConnectionsInterface* Connections,
                                   int _my_rank, int _num_ranks,
                                   unsigned int _receptor_type,
                                   unsigned int _external_function) {
  n_per_rank = (_neuron_num - 1) / _num_ranks + 1;
  int n_per_your_rank = (_neuron_num % _num_ranks == 0) ? n_per_rank:
                        (_my_rank <= _neuron_num/n_per_rank) ? n_per_rank:
                        (_my_rank == _neuron_num/n_per_rank + 1) ? _neuron_num % n_per_rank:
                        0;
  my_rank = _my_rank;
  num_ranks = _num_ranks;
  all_neuron_num = _neuron_num;

  first_neuron_id = n_per_rank * my_rank;
  last_neuron_id = n_per_rank * my_rank + n_per_your_rank - 1;
  neuron_num = last_neuron_id - first_neuron_id + 1;
  V_old = new double[neuron_num];
  V_new = new double[neuron_num];
  m = new double[neuron_num];
  h = new double[neuron_num];
  n = new double[neuron_num];
  I_ext = new double[neuron_num];
  I_syn = new double[neuron_num];

  recv_connections.resize(neuron_num);
  send_ids.resize(neuron_num);
  presynaptic_neurons = new remote_neuron*[neuron_num];

  Connections->GetAllNeuronIdsToSend(&send_ids);
  Connections->GetAllNeuronsToRecv(&recv_connections);

  communication_algorithm = (Connections->GetDensity() < collective_threshold) || (_neuron_num % _num_ranks != 0) ?
                            p2p : allgather;
#ifdef NEURODYNAMICS_DEBUG
  std::cout << "# Rank[" << my_rank << "] communication_algorithm=" << communication_algorithm <<
            " neuron_num=" << neuron_num <<
            " n_per_your_rank=" << n_per_your_rank << std::endl;
#endif
  MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_REMOTE_NEURON);
  MPI_Type_commit(&MPI_REMOTE_NEURON);

  if (communication_algorithm == p2p) {
    for (unsigned int i = 0; i < neuron_num; i++) {
      send_neuron_buffer = new remote_neuron[neuron_num];
      if (recv_connections[i].size() != 0)
        presynaptic_neurons[i] = new remote_neuron[recv_connections[i].size()];
      else
        presynaptic_neurons[i] = NULL;
    }
    recv_requests.resize(neuron_num);
    for (unsigned int iter_recv_neuron = 0; iter_recv_neuron < neuron_num; iter_recv_neuron++) {
      recv_requests[iter_recv_neuron].resize(recv_connections[iter_recv_neuron].size());
      for (unsigned int iter_neuron = 0; iter_neuron < recv_connections[iter_recv_neuron].size(); iter_neuron++) {
        MPI_Recv_init(&(presynaptic_neurons[iter_recv_neuron][iter_neuron]), 1, MPI_REMOTE_NEURON,
                      recv_connections[iter_recv_neuron][iter_neuron].id / n_per_rank,
                      recv_connections[iter_recv_neuron][iter_neuron].id,
                      MPI_COMM_WORLD, &(recv_requests[iter_recv_neuron][iter_neuron]));
      }
    }
    send_requests.resize(neuron_num);
    for (unsigned int iter_send_neuron = 0; iter_send_neuron < neuron_num; iter_send_neuron++) {
      send_requests[iter_send_neuron].resize(send_ids[iter_send_neuron].size());
      for (unsigned int iter_neuron = 0; iter_neuron < send_ids[iter_send_neuron].size(); iter_neuron++) {
        MPI_Send_init(&(send_neuron_buffer[iter_send_neuron]), 1, MPI_REMOTE_NEURON,
                      send_ids[iter_send_neuron][iter_neuron] / n_per_rank,
                      n_per_rank * my_rank + iter_send_neuron,
                      MPI_COMM_WORLD, &(send_requests[iter_send_neuron][iter_neuron]));
      }
    }
  } else {
    presynaptic_neurons[0] = new remote_neuron[all_neuron_num];
    for (unsigned int i = 1; i < neuron_num; i++) {
      presynaptic_neurons[i] = NULL;
    }
    send_neuron_buffer = NULL;
    allgather_request = new MPI_Request;
  }
  receptor_type = _receptor_type;
  /* External voltage init */
  switch (_external_function) {
    case I_EXTERNAL_NULL:
      for (unsigned int i = 0; i < neuron_num; i++) {
        I_ext[i] = 0.0;
      }
      break;

    case I_EXTERNAL_RANDOM:
      for (unsigned int i = 0; i < neuron_num; i++) {
        I_ext[i] = (static_cast<double>(rand()) / RAND_MAX)
                   * I_EXTERNAL_RANDOM_MAX_VALUE;
      }
      break;

    default:
      throw;
  }
  /* Starting variables init */
  for (unsigned int i = 0; i < neuron_num; i++) {
    V_old[i] = V_new[i] = V_REST;
    I_syn[i] = 0;

    m[i] = 0.055;
    h[i] = 0.59;
    n[i] = 0.32;
    if (communication_algorithm == p2p)
      send_neuron_buffer[i].Sact_generated = 0.0001;
  }
  if (my_rank == 0) {
    if (communication_algorithm == p2p)
    { std::cout << "Using point to point communications" << std::endl; }
    else
    { std::cout << "Using collective communications" << std::endl; }
  }
}

void NeuronHodgkinMPI::SendRecvPresynaptic() {
    if (send_requests[current_neuron].size() != 0)
    {
      send_neuron_buffer[current_neuron].Sact_generated =
              CalcPostSActivation(send_neuron_buffer[current_neuron].Sact_generated);
      send_neuron_buffer[current_neuron].voltage = V_old[current_neuron];
      MPI_Startall(send_requests[current_neuron].size(), send_requests[current_neuron].data());
    }
    if (recv_requests[current_neuron].size() != 0)
    {
      MPI_Startall(recv_requests[current_neuron].size(), recv_requests[current_neuron].data());
    }
}

int NeuronHodgkinMPI::process(int turns) {
  double* tmp_ptr_V;
  print(-1, OUTPUT_INIT, 3);
#ifdef NEURODYNAMICS_DEBUG
  std::cout << "Rank[" << my_rank << "] Starting process(" << turns << " iterations)." << std::endl;
#endif
  for (int cur_turn = 0; cur_turn < turns; cur_turn++) {
    if (communication_algorithm == p2p) {
      for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
        SendRecvPresynaptic();
      }
    }else {
      for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
        presynaptic_neurons[0][n_per_rank * my_rank + current_neuron].Sact_generated =
                CalcPostSActivation(presynaptic_neurons[0][n_per_rank * my_rank + current_neuron].Sact_generated);
        presynaptic_neurons[0][n_per_rank * my_rank + current_neuron].voltage = V_old[current_neuron];
      }
      MPI_Iallgather(MPI_IN_PLACE,  n_per_rank, MPI_REMOTE_NEURON,
                     &(presynaptic_neurons[0][0]), n_per_rank, MPI_REMOTE_NEURON,
                     MPI_COMM_WORLD, allgather_request);
      }
    for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
      m[current_neuron] = RungeKutta4(&NeuronHodgkinMPI::_CalcDerivativeFm, m[current_neuron], DELTAT);
      h[current_neuron] = RungeKutta4(&NeuronHodgkinMPI::_CalcDerivativeFh, h[current_neuron], DELTAT);
      n[current_neuron] = RungeKutta4(&NeuronHodgkinMPI::_CalcDerivativeFn, n[current_neuron], DELTAT);
    }
    if (communication_algorithm == p2p) {
      for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
        MPI_Waitall(send_requests[current_neuron].size(), send_requests[current_neuron].data(), MPI_STATUS_IGNORE);
      }
      for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
        MPI_Waitall(recv_requests[current_neuron].size(), recv_requests[current_neuron].data(), MPI_STATUS_IGNORE);
      }
      for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
        I_syn[current_neuron] = CalcSynapticCurrent();
      }
    } else {
      MPI_Wait(allgather_request, MPI_STATUS_IGNORE);
      for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
        I_syn[current_neuron] = CalcSynapticCurrentCollective();
      }
    }
    for (current_neuron = 0; current_neuron < neuron_num; current_neuron++) {
      V_new[current_neuron] = RungeKutta4(&NeuronHodgkinMPI::_CalcDerivativeFv, V_old[current_neuron], DELTAT);
    }
    tmp_ptr_V = V_new;
    V_new = V_old;
    V_old = tmp_ptr_V;
#if defined(OUTPUT_PRINT_STEP) && OUTPUT_PRINT_STEP >= 1
    if((cur_turn % OUTPUT_PRINT_STEP) == 0)
#endif
      print(cur_turn);
 }
#ifdef NEURODYNAMICS_DEBUG
  std::cout << "Rank[" << my_rank << "] Finishing process." << std::endl;
#endif
  print(-1, OUTPUT_FINAL, 3);
  return 0;
}

double NeuronHodgkinMPI::CalcSynapticCurrentCollective() {
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
              presynaptic_neurons[0][recv_connections[current_neuron][iter_neuron].id].Sact_generated *
              (V_old[current_neuron] - Esyn);
  }
  return result;
}


void NeuronHodgkinMPI::print(int step, std::string name, int process_level) {
  name.append("_");
  char my_rank_str[21]; // enough to hold all numbers up to 64-bits
  sprintf(my_rank_str, "%d", my_rank);
  name += my_rank_str;
  NeuronHodgkin::print(step, name, process_level);
}

NeuronHodgkinMPI::~NeuronHodgkinMPI() {
  if(communication_algorithm == allgather)
    delete allgather_request;
  else
    delete send_neuron_buffer;
}
