cmake_minimum_required(VERSION 3.2)
project(neurodynamics)

set(SOURCE_FILES
    connections_dense.cpp
    connections_dense_mpi.cpp
    main_neuro.cpp
    neuron_hodgin_mpi.cpp
    neuron_hodgin.cpp
    neuronetwork.cpp
)

add_executable(neurodynamics ${SOURCE_FILES})

if (MPI_CXX_LIBRARIES)
    target_link_libraries(neurodynamics ${MPI_CXX_LIBRARIES})
endif()

if(MPI_CXX_COMPILE_FLAGS)
  set_target_properties(neurodynamics PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_CXX_LINK_FLAGS)
  set_target_properties(neurodynamics PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()
