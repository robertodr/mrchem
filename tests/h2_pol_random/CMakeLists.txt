if(ENABLE_MPI)
  set(_h2_pol_random_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_static_polarizability_random"
  LABELS "H2_static_polarizability_random;polarizability;mrchem;h2_pol_lda"
  COST 100
  LAUNCH_AGENT ${_h2_pol_random_launcher}
  )
