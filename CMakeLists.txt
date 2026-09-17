cmake_minimum_required(VERSION 3.20)

project(FermiForge LANGUAGES Fortran)

include(CTest)

option(FERMIFORGE_ENABLE_MPI "Build the MPI field-map execution layer" ON)

add_library(he3_mesh
  src/he3_kinds.f90
  src/he3_quadrature.f90
  src/spin_matrix_2x2.f90
  src/quasiclassical_propagator.f90
  src/he3_self_consistency_integrand.f90
  src/order_parameter_basis.f90
  src/cartesian_mesh_2d.f90
  src/spinful_state_2d.f90
  src/spinful_field_io_2d.f90
  src/bilinear_field_sampler_2d.f90
  src/straight_trajectory_2d.f90
  src/trajectory_self_energy_2d.f90
  src/free_vortex_asymptotic_2d.f90
  src/riccati_segment.f90
  src/legacy_riccati_trajectory_2d.f90
  src/he3_trajectory_point_map.f90
  src/he3_serial_point_map.f90
  src/he3_field_map_2d.f90
  src/new_src_iteration_layout_2d.f90
  src/legacy_anderson_mixing.f90
  src/he3_nonlinear_iteration_2d.f90
  src/radial_mesh.f90
  src/radial_interpolation.f90
  src/radial_state.f90
  src/axial_radial_embedding.f90
  src/legacy_radial_profile_io.f90
)

set_target_properties(he3_mesh PROPERTIES
  Fortran_STANDARD 2018
  Fortran_STANDARD_REQUIRED YES
  Fortran_EXTENSIONS NO
)

add_executable(fermiforge_2d_demo app/fermiforge_2d_demo.f90)
target_link_libraries(fermiforge_2d_demo PRIVATE he3_mesh)
set_target_properties(fermiforge_2d_demo PROPERTIES
  Fortran_STANDARD 2018
  Fortran_STANDARD_REQUIRED YES
  Fortran_EXTENSIONS NO
)

add_executable(compare_radial_point_map app/compare_radial_point_map.f90)
target_link_libraries(compare_radial_point_map PRIVATE he3_mesh)
set_target_properties(compare_radial_point_map PROPERTIES
  Fortran_STANDARD 2018
  Fortran_STANDARD_REQUIRED YES
  Fortran_EXTENSIONS NO
)

if(FERMIFORGE_ENABLE_MPI)
  find_package(MPI COMPONENTS Fortran QUIET)
  if(MPI_Fortran_FOUND)
    add_library(he3_mpi_execution src/he3_mpi_field_map_2d.f90)
    target_link_libraries(he3_mpi_execution PUBLIC he3_mesh MPI::MPI_Fortran)
    set_target_properties(he3_mpi_execution PROPERTIES
      Fortran_STANDARD 2018
      Fortran_STANDARD_REQUIRED YES
      Fortran_EXTENSIONS NO
    )
    add_executable(fermiforge_2d_mpi_smoke
      app/fermiforge_2d_mpi_smoke.f90
    )
    target_link_libraries(fermiforge_2d_mpi_smoke PRIVATE he3_mpi_execution)
    set_target_properties(fermiforge_2d_mpi_smoke PROPERTIES
      Fortran_STANDARD 2018
      Fortran_STANDARD_REQUIRED YES
      Fortran_EXTENSIONS NO
    )
    add_executable(benchmark_axisymmetric_core_2d
      app/benchmark_axisymmetric_core_2d.f90
    )
    target_link_libraries(benchmark_axisymmetric_core_2d
      PRIVATE he3_mpi_execution
    )
    set_target_properties(benchmark_axisymmetric_core_2d PROPERTIES
      Fortran_STANDARD 2018
      Fortran_STANDARD_REQUIRED YES
      Fortran_EXTENSIONS NO
    )
    set(FERMIFORGE_MPI_1_LOCAL_FLAGS)
    set(FERMIFORGE_MPI_4_LOCAL_FLAGS)
    set(FERMIFORGE_MPI_10_LOCAL_FLAGS)
    if(APPLE)
      execute_process(
        COMMAND ${MPIEXEC_EXECUTABLE} --version
        OUTPUT_VARIABLE FERMIFORGE_MPI_VERSION_TEXT
        ERROR_VARIABLE FERMIFORGE_MPI_VERSION_TEXT
      )
      if(FERMIFORGE_MPI_VERSION_TEXT MATCHES "Open MPI")
        # Open MPI 5 may be unable to infer macOS core slots in restricted
        # test environments. These local smoke tests declare their slots
        # explicitly and do not affect production launch commands.
        set(FERMIFORGE_MPI_1_LOCAL_FLAGS
            --host localhost:1 --map-by ppr:1:node --bind-to none)
        set(FERMIFORGE_MPI_4_LOCAL_FLAGS
            --host localhost:4 --map-by ppr:4:node --bind-to none)
        set(FERMIFORGE_MPI_10_LOCAL_FLAGS
            --host localhost:10 --map-by ppr:10:node --bind-to none)
      endif()
    endif()
  else()
    message(STATUS "MPI Fortran not found; MPI field-map layer is disabled")
  endif()
endif()

if(BUILD_TESTING)
  add_executable(test_radial_mesh tests/test_radial_mesh.f90)
  target_link_libraries(test_radial_mesh PRIVATE he3_mesh)
  set_target_properties(test_radial_mesh PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME radial_mesh COMMAND test_radial_mesh)

  add_executable(test_radial_state tests/test_radial_state.f90)
  target_link_libraries(test_radial_state PRIVATE he3_mesh)
  set_target_properties(test_radial_state PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME radial_state COMMAND test_radial_state)

  add_executable(test_legacy_anderson_mixing
    tests/test_legacy_anderson_mixing.f90
  )
  target_link_libraries(test_legacy_anderson_mixing PRIVATE he3_mesh)
  set_target_properties(test_legacy_anderson_mixing PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME legacy_anderson_mixing COMMAND test_legacy_anderson_mixing)

  add_executable(test_cartesian_field_2d
    tests/test_cartesian_field_2d.f90
  )
  target_link_libraries(test_cartesian_field_2d PRIVATE he3_mesh)
  set_target_properties(test_cartesian_field_2d PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME cartesian_field_2d COMMAND test_cartesian_field_2d)

  add_executable(test_order_parameter_basis
    tests/test_order_parameter_basis.f90
  )
  target_link_libraries(test_order_parameter_basis PRIVATE he3_mesh)
  set_target_properties(test_order_parameter_basis PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME order_parameter_basis COMMAND test_order_parameter_basis)

  add_executable(test_axial_radial_embedding
    tests/test_axial_radial_embedding.f90
  )
  target_link_libraries(test_axial_radial_embedding PRIVATE he3_mesh)
  set_target_properties(test_axial_radial_embedding PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME axial_radial_embedding COMMAND test_axial_radial_embedding)

  add_executable(test_straight_trajectory_2d
    tests/test_straight_trajectory_2d.f90
  )
  target_link_libraries(test_straight_trajectory_2d PRIVATE he3_mesh)
  set_target_properties(test_straight_trajectory_2d PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME straight_trajectory_2d COMMAND test_straight_trajectory_2d)

  add_executable(test_riccati_segment
    tests/test_riccati_segment.f90
  )
  target_link_libraries(test_riccati_segment PRIVATE he3_mesh)
  set_target_properties(test_riccati_segment PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME riccati_segment COMMAND test_riccati_segment)

  add_executable(test_trajectory_self_energy_2d
    tests/test_trajectory_self_energy_2d.f90
  )
  target_link_libraries(test_trajectory_self_energy_2d PRIVATE he3_mesh)
  set_target_properties(test_trajectory_self_energy_2d PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME trajectory_self_energy_2d COMMAND test_trajectory_self_energy_2d)

  add_executable(test_free_vortex_asymptotic_2d
    tests/test_free_vortex_asymptotic_2d.f90
  )
  target_link_libraries(test_free_vortex_asymptotic_2d PRIVATE he3_mesh)
  set_target_properties(test_free_vortex_asymptotic_2d PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME free_vortex_asymptotic_2d
           COMMAND test_free_vortex_asymptotic_2d)

  add_executable(test_legacy_riccati_trajectory_2d
    tests/test_legacy_riccati_trajectory_2d.f90
  )
  target_link_libraries(test_legacy_riccati_trajectory_2d PRIVATE he3_mesh)
  set_target_properties(test_legacy_riccati_trajectory_2d PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME legacy_riccati_trajectory_2d
           COMMAND test_legacy_riccati_trajectory_2d)

  add_executable(test_spin_matrix_propagator
    tests/test_spin_matrix_propagator.f90
  )
  target_link_libraries(test_spin_matrix_propagator PRIVATE he3_mesh)
  set_target_properties(test_spin_matrix_propagator PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME spin_matrix_propagator COMMAND test_spin_matrix_propagator)

  add_executable(test_he3_self_consistency_integrand
    tests/test_he3_self_consistency_integrand.f90
  )
  target_link_libraries(test_he3_self_consistency_integrand PRIVATE he3_mesh)
  set_target_properties(test_he3_self_consistency_integrand PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME he3_self_consistency_integrand
           COMMAND test_he3_self_consistency_integrand)

  add_executable(test_he3_trajectory_point_map
    tests/test_he3_trajectory_point_map.f90
  )
  target_link_libraries(test_he3_trajectory_point_map PRIVATE he3_mesh)
  set_target_properties(test_he3_trajectory_point_map PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME he3_trajectory_point_map COMMAND test_he3_trajectory_point_map)

  add_executable(test_he3_quadrature tests/test_he3_quadrature.f90)
  target_link_libraries(test_he3_quadrature PRIVATE he3_mesh)
  set_target_properties(test_he3_quadrature PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(
    NAME he3_quadrature
    COMMAND test_he3_quadrature
            ${CMAKE_CURRENT_SOURCE_DIR}/new_src/gauss11.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/new_src/ozaki.dat
  )

  add_executable(test_he3_serial_point_map tests/test_he3_serial_point_map.f90)
  target_link_libraries(test_he3_serial_point_map PRIVATE he3_mesh)
  set_target_properties(test_he3_serial_point_map PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(
    NAME he3_serial_point_map
    COMMAND test_he3_serial_point_map
            ${CMAKE_CURRENT_SOURCE_DIR}/new_src/gauss11.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/new_src/ozaki.dat
  )

  add_executable(test_legacy_radial_profile_io
    tests/test_legacy_radial_profile_io.f90
  )
  target_link_libraries(test_legacy_radial_profile_io PRIVATE he3_mesh)
  set_target_properties(test_legacy_radial_profile_io PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(
    NAME legacy_radial_profile_io
    COMMAND test_legacy_radial_profile_io
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/a_phase_core_2d/reference/current_new_src_T0.30_F1s5.4/op_xyz
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/a_phase_core_2d/reference/current_new_src_T0.30_F1s5.4/curr
            ${CMAKE_CURRENT_SOURCE_DIR}/new_src/op_xyz_n
            ${CMAKE_CURRENT_SOURCE_DIR}/new_src/curr_n
  )

  add_executable(test_he3_field_map_2d tests/test_he3_field_map_2d.f90)
  target_link_libraries(test_he3_field_map_2d PRIVATE he3_mesh)
  set_target_properties(test_he3_field_map_2d PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
  )
  add_test(NAME he3_field_map_2d COMMAND test_he3_field_map_2d)

  if(TARGET he3_mpi_execution)
    add_executable(test_he3_mpi_field_map tests/test_he3_mpi_field_map.f90)
    target_link_libraries(test_he3_mpi_field_map PRIVATE he3_mpi_execution)
    set_target_properties(test_he3_mpi_field_map PROPERTIES
      Fortran_STANDARD 2018
      Fortran_STANDARD_REQUIRED YES
      Fortran_EXTENSIONS NO
    )
    add_test(
      NAME he3_mpi_field_map_1_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${FERMIFORGE_MPI_1_LOCAL_FLAGS}
              ${MPIEXEC_NUMPROC_FLAG} 1
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_he3_mpi_field_map>
              ${CMAKE_CURRENT_SOURCE_DIR}/new_src/gauss11.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/new_src/ozaki.dat
              ${MPIEXEC_POSTFLAGS}
    )
    add_test(
      NAME he3_mpi_field_map_4_ranks
      COMMAND ${MPIEXEC_EXECUTABLE} ${FERMIFORGE_MPI_4_LOCAL_FLAGS}
              ${MPIEXEC_NUMPROC_FLAG} 4
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_he3_mpi_field_map>
              ${CMAKE_CURRENT_SOURCE_DIR}/new_src/gauss11.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/new_src/ozaki.dat
              ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      he3_mpi_field_map_1_rank PROPERTIES PROCESSORS 1
    )
    set_tests_properties(
      he3_mpi_field_map_4_ranks PROPERTIES PROCESSORS 4
    )
    add_test(
      NAME he3_mpi_field_map_10_ranks
      COMMAND ${MPIEXEC_EXECUTABLE} ${FERMIFORGE_MPI_10_LOCAL_FLAGS}
              ${MPIEXEC_NUMPROC_FLAG} 10
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_he3_mpi_field_map>
              ${CMAKE_CURRENT_SOURCE_DIR}/new_src/gauss11.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/new_src/ozaki.dat
              ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      he3_mpi_field_map_10_ranks PROPERTIES PROCESSORS 10
    )
  endif()

  add_test(
    NAME a_phase_core_radial_point_map_comparison
    COMMAND compare_radial_point_map
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/a_phase_core_2d/reference/current_new_src_T0.30_F1s5.4/op_xyz
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/a_phase_core_2d/reference/current_new_src_T0.30_F1s5.4/curr
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/a_phase_core_2d/reference/current_new_src_T0.30_F1s5.4/gauss11.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/a_phase_core_2d/reference/current_new_src_T0.30_F1s5.4/ozaki.dat
            512 32
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/a_phase_core_2d/reference/current_new_src_T0.30_F1s5.4/point_24.expected
  )

  add_test(
    NAME normal_core_radial_point_map_comparison
    COMMAND compare_radial_point_map
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/normal_core_2d/reference/current_new_src_T0.30_F1s5.4/op_xyz
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/normal_core_2d/reference/current_new_src_T0.30_F1s5.4/curr
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/normal_core_2d/reference/current_new_src_T0.30_F1s5.4/gauss11.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/normal_core_2d/reference/current_new_src_T0.30_F1s5.4/ozaki.dat
            512 32
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/normal_core_2d/reference/current_new_src_T0.30_F1s5.4/point_24.expected
  )

  add_test(
    NAME fermiforge_2d_demo
    COMMAND fermiforge_2d_demo
            ${CMAKE_CURRENT_SOURCE_DIR}/examples/2d_sampler_demo.nml
            ${CMAKE_CURRENT_BINARY_DIR}/fermiforge_2d_demo_test.dat
  )
endif()
