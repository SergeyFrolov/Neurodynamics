// Copyright 2015 Sergey Frolov. All rights reserved.
// Use of this source code is governed by a LGPL license that can be
// found in the LICENSE file.

#ifndef NEURODYNAMICS_DEFINES_H_
#define NEURODYNAMICS_DEFINES_H_

#define ind(x, y, width) ((x) + (y) * (width))
#define MAX(x, y) (x) > (y) ? (x) : (y)

/* run config */
#define NEURON_NUM     100
#define PROCESS_TURNS  200000
#define RANDOMIZE_SEED false

#ifndef MPI_COLLECTIVE_DENSITY_THRESHOLD
#define MPI_COLLECTIVE_DENSITY_THRESHOLD 0.7
// 0 - empty graph, 1 - full
// Graphs with density higher than threshold would be processed by collective algorithms
#endif

/* topology-related */
#define GRAPH_RAND                  0x00000010
#define GRAPH_RAND_P_CONNECTION     0.00000000     // probability of connection in random
#define GRAPH_RING                  0x00000020
#define GRAPH_CONNECT_EPS           0.000001
//  connections weaker than CONNECT_EPS would be disregarded
#define GRAPH_DEFAULT               GRAPH_RAND

#define DELTAT                      0.001  //  step

/* external current */
#define I_EXTERNAL_NULL             0x00000000
#define I_EXTERNAL_RANDOM           0x00000020
#define I_EXTERNAL_UNIFORM          0x00000040
#define I_EXTERNAL_DEFAULT          I_EXTERNAL_UNIFORM

#define I_EXTERNAL_MIN_VALUE 6.3
#define I_EXTERNAL_MAX_VALUE 18.
// values would be in following range:
// [0, I_EXTERNAL_MAX_VALUE]
// could be negative

/* receptor-related */
#define AMPA_RECEPTOR     0x00000010
#define NMDA_RECEPTOR     0x00000020
#define GABA_A_RECEPTOR   0x00000030
#define GABA_B_RECEPTOR   0x00000031
#define DEFAULT_RECEPTOR  AMPA_RECEPTOR

#define gAMPA_MIN         0.35
#define gAMPA_MAX         1.00
#define gNMDA_MIN         0.01
#define gNMDA_MAX         0.60
#define gGABAA_MIN        0.25
#define gGABAA_MAX        1.20
#define gGABAB_MIN        0.90
#define gGABAB_MAX        1.10

#define AMPA_E            0.       // mV
// reversal potential
#define AMPA_ALPHA        1100     // mM^(-1)s^(-1)
#define AMPA_BETA         190      //        s^(-1)
// alpha and beta are voltage-independent forward and backward rate constants

#define RECEPTOR_Kp       5.     // mV
// steepness
#define RECEPTOR_Vp       2.     // mV
// the value at which the function is half - activated.
#define RECEPTOR_Tmax     1.     // mM
// max concentration of transmitter in the synaptic cleft

/* model-related */
#define C_M               1.    // uF/cmI

#define V_REST            -75.                  // mV
#define V_NA_REST         (V_REST + 115)        // mV
#define V_K_REST          (V_REST - 12)         // mV
#define V_L_REST          (V_REST + 10.6130)    // mV

#define G_NA_MAX          120  // ms/cmI
#define G_K_MAX           36   // ms/cmI
#define G_L               0.3  // ms/cmI

/* output-related */
#define OUTPUT_TOPOLOGY        "topology"
#define OUTPUT_INIT            "init_state"  // including I_ext
#define OUTPUT_FINAL           "final_state"
#define OUTPUT_PROCESS         "mid_state"
#define OUTPUT_SKIP_TURNS      100000
#ifndef OUTPUT_PROCESS_LEVEL
#define OUTPUT_PROCESS_LEVEL   0
#endif
#ifndef OUTPUT_PRINT_STEP
#define OUTPUT_PRINT_STEP      1 // 1 - every step, 10 - every 10th
#endif
// 0 - no mid output (only init and final)
// 1 - only peaks
// 2 - only V_new
// 3 - everything, but Sact
// 4 - everything and Sact

#define NEURON_PEAK_REGISTER_THRESHOLD 0  // mV
#define NEURON_PEAK_COOLDOWN_THRESHOLD -10  // mV
// neuron should have V > REGISTER_THRESHOLD to register the peak
// after it, neuron should get voltage smaller than COOLDOWN_THRESHOLD to register new peak

// external defines:
// NEURODYNAMICS_WITH_MPI
// NEURODYNAMICS_DEBUG
// OUTPUT_VERBOSE

#endif  // NEURODYNAMICS_DEFINES_H_
