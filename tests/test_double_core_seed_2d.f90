program test_double_core_seed_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
    make_uniform_cartesian_mesh
  use double_core_seed_2d, only : double_core_seed_report_t, &
    initialize_regularized_london_double_core_2d
  use spinful_state_2d, only : spinful_state_2d_t
  implicit none

  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state
  type(double_core_seed_report_t) :: report
  complex(rk) :: center_xz, center_zx
  integer :: center, positive_core, negative_core
  real(rk), parameter :: gap = 0.28_rk
  real(rk), parameter :: half_offset = 2.0_rk
  real(rk), parameter :: tolerance = 2.0e-13_rk

  call make_uniform_cartesian_mesh(-12.0_rk, 12.0_rk, 48, &
                                   -12.0_rk, 12.0_rk, 48, mesh)
  call initialize_regularized_london_double_core_2d( &
    mesh, gap, half_offset, 0.8_rk, 0.6_rk, state, report)

  call require(state%is_valid_for(mesh), "London seed state is invalid")
  call require(maxval(abs(state%current_mean_field)) <= tiny(1.0_rk), &
    "London seed mean field is not zero")

  center = mesh%point_index(24, 24)
  positive_core = mesh%point_index(24, 28)
  negative_core = mesh%point_index(24, 20)
  center_xz = state%order_parameter(1, 3, center)
  center_zx = state%order_parameter(3, 1, center)
  call require(abs(center_xz + center_zx) < tolerance, &
    "London seed lost its antisymmetric xz/zx planar component")
  call require(abs(center_xz) > 0.1_rk * gap, &
    "London seed has no planar core filling")
  call require(abs(state%order_parameter(2, 2, center)) < tolerance, &
    "London seed domain wall does not suppress A_yy at the center")

  call require(abs(state%order_parameter(2, 2, positive_core)) < tolerance, &
    "positive half core is not on the regularized domain wall")
  call require(abs(state%order_parameter(2, 2, negative_core)) < tolerance, &
    "negative half core is not on the regularized domain wall")
  call require(report%positive_half_core_pair_amplitude < &
               report%center_pair_amplitude, &
    "positive half core is not a pair-amplitude minimum")
  call require(report%negative_half_core_pair_amplitude < &
               report%center_pair_amplitude, &
    "negative half core is not a pair-amplitude minimum")
  call require(report%outer_pair_amplitude > report%center_pair_amplitude, &
    "regularized London seed does not recover toward the bulk")
  call require(abs(report%positive_half_core_pair_amplitude - &
                   report%negative_half_core_pair_amplitude) < tolerance, &
    "regularized London seed does not place symmetric half cores")

  print '(a)', "double-core seed tests passed"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_double_core_seed_2d
