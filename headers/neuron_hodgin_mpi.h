// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#ifndef NEURODYNAMICS_NEURON_HODGIN_MPI_H_
#define NEURODYNAMICS_NEURON_HODGIN_MPI_H_

#pragma once

#include "neuron_hodgin.h"
#include "mpi.h"

class NeuronHodgkinMPI : public NeuronHodgkin {
  enum algorithm {p2p, allgather};
  double collective_threshold = MPI_COLLECTIVE_DENSITY_THRESHOLD;
  algorithm communication_algorithm;

  int my_rank;
  int num_ranks;
  int all_neuron_num;
  int n_per_rank;

  int first_neuron_id;
  int last_neuron_id;

  MPI_Datatype MPI_REMOTE_NEURON;
  remote_neuron* send_neuron_buffer;

  std::vector<std::vector<MPI_Request>> recv_requests;
  std::vector<std::vector<MPI_Request>> send_requests;
  MPI_Request* allgather_request;
 public:
  NeuronHodgkinMPI(unsigned int _neuron_num, ConnectionsInterface* Connections,
                   int _my_rank, int _num_ranks,
                   unsigned int _receptor_type = AMPA_RECEPTOR,
                   unsigned int _external_function = I_EXTERNAL_NULL);
  int process(int turns) override;
  void SendRecvPresynaptic() override;
  void print(int step, std::string name = OUTPUT_PROCESS, int process_level = OUTPUT_PROCESS_LEVEL) override;
  double CalcSynapticCurrentCollective();
  ~NeuronHodgkinMPI();
};
#endif  // NEURODYNAMICS_NEURON_HODGIN_MPI_H_
