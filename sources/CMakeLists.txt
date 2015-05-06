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

add_executable(${TARGET_NAME} ${SOURCE_FILES})

if (MPI_CXX_LIBRARIES)
    target_link_libraries(${TARGET_NAME} ${MPI_CXX_LIBRARIES})
endif()

if(MPI_CXX_COMPILE_FLAGS)
  set_target_properties(${TARGET_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_CXX_LINK_FLAGS)
  set_target_properties(${TARGET_NAME} PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()
