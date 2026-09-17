program compare_radial_point_map
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t
  use axial_radial_embedding, only : axial_radial_profile_t, &
                                     embed_axial_radial_profile_2d
  use legacy_radial_profile_io, only : read_new_src_radial_profile
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
                                make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t, &
                             read_legacy_gauss_table, &
                             read_legacy_ozaki_table
  use he3_serial_point_map, only : he3_point_map_diagnostics_t, &
                                   evaluate_serial_he3_point_map
  implicit none

  real(rk), parameter :: source_feedback = 5.4_rk / (1.0_rk + 5.4_rk / 3.0_rk)
  real(rk), parameter :: source_step = 0.175_rk
  real(rk), parameter :: source_boundary_relaxation = 10.0_rk * source_step

  type(angular_quadrature_3d_t) :: angular
  type(axial_radial_profile_t) :: profile
  type(cartesian_mesh_2d_t) :: cartesian_grid
  type(he3_point_map_diagnostics_t) :: diagnostics
  type(ozaki_quadrature_t) :: energy
  type(radial_mesh_t) :: radial_grid
  type(spinful_state_2d_t) :: state
  complex(rk) :: mapped_gap(3, 3), reference_gap(3, 3)
  real(rk) :: mapped_current(3), reference_current(3)
  real(rk) :: gap_absolute_difference, gap_relative_difference, half_width
  real(rk) :: mean_field_absolute_difference, reference_radius
  character(len=1024) :: argument, current_path, gauss_path, order_parameter_path
  character(len=1024) :: ozaki_path, reference_path
  integer :: cartesian_cells, component, orbital, radial_index, spin
  logical :: enforce_reference_tolerance, succeeded

  if (command_argument_count() < 4 .or. command_argument_count() > 7) &
    error stop "usage: compare_radial_point_map op_xyz curr gauss11.dat " // &
               "ozaki.dat [cartesian_cells [half_width [expected_file]]]"
  call get_command_argument(1, order_parameter_path)
  call get_command_argument(2, current_path)
  call get_command_argument(3, gauss_path)
  call get_command_argument(4, ozaki_path)
  cartesian_cells = 256
  half_width = 32.0_rk
  radial_index = 24
  reference_path = ""
  if (command_argument_count() >= 5) then
    call get_command_argument(5, argument)
    read(argument, *) cartesian_cells
  end if
  if (command_argument_count() >= 6) then
    call get_command_argument(6, argument)
    read(argument, *) half_width
  end if
  if (command_argument_count() >= 7) &
    call get_command_argument(7, reference_path)
  if (cartesian_cells < 4 .or. half_width <= 0.0_rk) &
    error stop "invalid Cartesian comparison mesh"

  call read_new_src_radial_profile( &
    trim(order_parameter_path), trim(current_path), radial_grid, profile)
  call read_legacy_gauss_table(trim(gauss_path), 48, 11, angular)
  call read_legacy_ozaki_table(trim(ozaki_path), energy)
  call make_uniform_cartesian_mesh( &
    -half_width, half_width, cartesian_cells, &
    -half_width, half_width, cartesian_cells, cartesian_grid)
  call allocate_spinful_state_2d(cartesian_grid, state)
  call embed_axial_radial_profile_2d( &
    radial_grid, profile, 1.0_rk, cartesian_grid, state)

  reference_gap = reshape([ &
    cmplx( 6.796686615524666e-2_rk, -1.833571904020656e-12_rk, rk), &
    cmplx( 9.854057500716621e-13_rk,  2.200015455446118e-2_rk, rk), &
    cmplx( 3.010035799017428e-1_rk, -7.924273366690986e-12_rk, rk), &
    cmplx(-8.768469192899934e-13_rk, -1.831575239097041e-2_rk, rk), &
    cmplx( 4.405380514573941e-2_rk, -9.495215364795768e-13_rk, rk), &
    cmplx( 6.635792903504480e-12_rk,  2.497188329805768e-1_rk, rk), &
    cmplx(-1.704971548761743e-1_rk,  4.478350116529464e-12_rk, rk), &
    cmplx(-8.262891126590807e-13_rk, -2.369278046987337e-2_rk, rk), &
    cmplx( 8.803873255759169e-2_rk, -2.588614107993685e-12_rk, rk)], &
    [3, 3])
  reference_current = [ &
     8.826461765941074e-16_rk, -1.679384493624108e-2_rk, &
     2.115086769252801e-14_rk]
  reference_radius = 3.6095394702719363_rk
  if (len_trim(reference_path) > 0) &
    call read_reference_file( &
      trim(reference_path), reference_gap, reference_current, radial_index, &
      reference_radius)
  if (radial_index < 0 .or. radial_index >= radial_grid%point_count()) &
    error stop "reference radial index lies outside the source profile"
  if (abs(radial_grid%r(radial_index) - reference_radius) > 1.0e-12_rk) &
    error stop "reference radius does not match the reconstructed radial mesh"

  call evaluate_serial_he3_point_map( &
    cartesian_grid, state, [radial_grid%r(radial_index), 0.0_rk], &
    angular, energy, source_feedback, source_step, 1, &
    source_boundary_relaxation, mapped_gap, mapped_current, diagnostics, &
    succeeded)
  if (.not. succeeded) error stop "modern radial point-map evaluation failed"

  gap_absolute_difference = maxval(abs(mapped_gap - reference_gap))
  gap_relative_difference = gap_absolute_difference / maxval(abs(reference_gap))
  mean_field_absolute_difference = &
    maxval(abs(mapped_current - reference_current))

  print '(a,l1)', "succeeded: ", succeeded
  print '(a,i0)', "directions: ", diagnostics%direction_count
  print '(a,i0)', "poles: ", diagnostics%pole_count
  print '(a,i0)', "Cartesian cells per axis: ", cartesian_cells
  print '(a,es12.4)', "Cartesian spacing: ", &
    2.0_rk * half_width / real(cartesian_cells, rk)
  print '(a,i0)', "accepted contributions: ", &
    diagnostics%accepted_contributions
  print '(a,es12.4)', "maximum normalization error: ", &
    diagnostics%maximum_normalization_error
  print '(a,3(1x,f8.3))', "timing [sample, propagate, accumulate] s:", &
    diagnostics%sampling_seconds, diagnostics%propagation_seconds, &
    diagnostics%accumulation_seconds
  print '(a,es12.4)', "maximum absolute gap difference: ", &
    gap_absolute_difference
  print '(a,es12.4)', "maximum relative gap difference: ", &
    gap_relative_difference
  print '(a,es12.4)', "maximum absolute mean-field difference: ", &
    mean_field_absolute_difference
  do spin = 1, 3
    do orbital = 1, 3
      component = (spin - 1) * 3 + orbital
      print '(a,1x,i2,4(1x,es18.10))', "MAP_COMPONENT", component, &
        real(mapped_gap(spin, orbital), rk), &
        aimag(mapped_gap(spin, orbital)), &
        real(reference_gap(spin, orbital), rk), &
        aimag(reference_gap(spin, orbital))
    end do
  end do
  do component = 1, 3
    print '(a,1x,i2,2(1x,es18.10))', "MAP_MEAN_FIELD", component, &
      mapped_current(component), reference_current(component)
  end do
  enforce_reference_tolerance = command_argument_count() == 4 .or. &
    len_trim(reference_path) > 0
  if (enforce_reference_tolerance) then
    if (gap_relative_difference >= 1.0e-3_rk) &
      error stop "default radial point-map gap tolerance was not met"
    if (mean_field_absolute_difference >= 5.0e-5_rk) &
      error stop "default radial point-map mean-field tolerance was not met"
  end if

contains

  subroutine read_reference_file( &
      path, gap, current, reference_point, radius)
    character(len=*), intent(in) :: path
    complex(rk), intent(out) :: gap(3, 3)
    real(rk), intent(out) :: current(3)
    integer, intent(out) :: reference_point
    real(rk), intent(out) :: radius

    character(len=4096) :: line
    character(len=64) :: label
    logical :: component_seen(12), index_seen, radius_seen
    real(rk) :: imaginary_part, real_part
    integer :: component_index, input_status, unit

    gap = cmplx(0.0_rk, 0.0_rk, rk)
    current = 0.0_rk
    component_seen = .false.
    index_seen = .false.
    radius_seen = .false.
    open(newunit=unit, file=trim(path), status="old", action="read", &
         iostat=input_status)
    if (input_status /= 0) error stop "could not open point-map reference file"
    do
      read(unit, '(a)', iostat=input_status) line
      if (input_status < 0) exit
      if (input_status > 0) error stop "could not read point-map reference file"
      if (index(line, "REFERENCE_RADIAL_INDEX") == 1) then
        read(line, *, iostat=input_status) label, reference_point
        if (input_status /= 0) error stop "invalid reference radial index"
        index_seen = .true.
      else if (index(line, "REFERENCE_RADIUS") == 1) then
        read(line, *, iostat=input_status) label, radius
        if (input_status /= 0) error stop "invalid reference radius"
        radius_seen = .true.
      else if (index(line, "REFERENCE_COMPONENT") == 1) then
        read(line, *, iostat=input_status) label, component_index, &
          real_part, imaginary_part
        if (input_status /= 0 .or. component_index < 1 .or. &
            component_index > 12) &
          error stop "invalid point-map reference component"
        if (component_seen(component_index)) &
          error stop "duplicate point-map reference component"
        component_seen(component_index) = .true.
        if (component_index <= 9) then
          gap((component_index - 1) / 3 + 1, &
              modulo(component_index - 1, 3) + 1) = &
            cmplx(real_part, imaginary_part, rk)
        else
          current(component_index - 9) = real_part
          if (abs(imaginary_part) > 1.0e-14_rk) &
            error stop "reference current component is unexpectedly complex"
        end if
      end if
    end do
    close(unit)
    if (.not. index_seen .or. .not. radius_seen .or. &
        .not. all(component_seen)) &
      error stop "point-map reference file is incomplete"
  end subroutine read_reference_file
end program compare_radial_point_map
