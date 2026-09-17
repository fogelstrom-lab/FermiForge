program test_legacy_radial_profile_io
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t
  use axial_radial_embedding, only : axial_radial_profile_t, &
                                     sample_axial_radial_profile
  use legacy_radial_profile_io, only : read_new_src_radial_profile
  implicit none

  type(radial_mesh_t) :: mesh
  type(axial_radial_profile_t) :: profile
  complex(rk) :: order_parameter(3, 3)
  real(rk) :: current_mean_field(3)
  character(len=16) :: grid_kind
  character(len=1024) :: current_path, normal_current_path
  character(len=1024) :: normal_order_parameter_path, order_parameter_path
  logical :: inside

  if (command_argument_count() /= 4) &
    error stop "radial profile I/O test needs tangent and uniform state paths"
  call get_command_argument(1, order_parameter_path)
  call get_command_argument(2, current_path)
  call get_command_argument(3, normal_order_parameter_path)
  call get_command_argument(4, normal_current_path)
  call read_new_src_radial_profile( &
    trim(order_parameter_path), trim(current_path), mesh, profile, &
    detected_grid_kind=grid_kind)

  call require(trim(grid_kind) == "tangent", &
               "current new_src tangent grid was not detected")
  call require(mesh%point_count() == 100, "new_src radial point count changed")
  call require(abs(mesh%r(24) - 3.6095394702719363_rk) < 6.0e-4_rk, &
               "new_src radial coordinate was not read correctly")
  call sample_axial_radial_profile( &
    mesh, profile, mesh%r(24), 0.0_rk, 1.0_rk, order_parameter, &
    current_mean_field, inside)
  call require(inside, "loaded radial node was rejected")
  call require(abs(current_mean_field(1)) < 1.0e-14_rk .and. &
               abs(current_mean_field(2) + 1.67938448e-2_rk) < 2.0e-9_rk .and. &
               abs(current_mean_field(3)) < 1.0e-14_rk, &
               "new_src azimuthal mean field was not loaded correctly")
  call require(maxval(abs(order_parameter)) > 0.1_rk, &
               "loaded new_src order parameter is unexpectedly empty")

  call read_new_src_radial_profile( &
    trim(normal_order_parameter_path), trim(normal_current_path), mesh, &
    profile, radial_grid_kind="auto", detected_grid_kind=grid_kind)
  call require(trim(grid_kind) == "uniform", &
               "archived normal-core uniform grid was not detected")
  call require(mesh%point_count() == 50, &
               "archived normal-core point count changed")
  call require(abs(mesh%r(24) - 30.0_rk * 24.0_rk / 49.0_rk) < 1.0e-13_rk, &
               "normal-core uniform coordinate was not reconstructed")
  call sample_axial_radial_profile( &
    mesh, profile, mesh%r(24), 0.0_rk, 1.0_rk, order_parameter, &
    current_mean_field, inside)
  call require(inside .and. maxval(abs(order_parameter)) > 0.1_rk, &
               "archived normal-core profile was not loaded")

  print '(a)', "new_src tangent and archived uniform profile I/O tests passed"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_legacy_radial_profile_io
