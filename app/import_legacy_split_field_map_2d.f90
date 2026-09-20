program import_legacy_split_field_map_2d
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t
  use spinful_field_io_2d, only : write_spinful_field_map_2d
  use legacy_split_field_io_2d, only : &
    legacy_split_field_report_2d_t, &
    read_legacy_split_field_directory_2d
  implicit none

  character(len=2048) :: directory, output_path
  integer :: argument_count
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state
  type(legacy_split_field_report_2d_t) :: report

  argument_count = command_argument_count()
  if (argument_count < 1 .or. argument_count > 2) then
    error stop "usage: import_legacy_split_field_map_2d DIRECTORY [OUTPUT]"
  end if
  call get_command_argument(1, directory)
  call read_legacy_split_field_directory_2d( &
    trim(directory), mesh, state, report)

  print '(a,i0,a,i0,a)', "grid: ", report%number_of_x_points, " x ", &
    report%number_of_y_points, " points"
  print '(a,i0)', "total points: ", report%number_of_points
  print '(a,es14.6)', "maximum stored gap-norm difference: ", &
    report%maximum_gap_norm_difference
  print '(a,es14.6)', "maximum stored in-plane magnitude difference: ", &
    report%maximum_in_plane_magnitude_difference
  print '(a,es14.6,a,2(es14.6,1x))', "minimum gap norm: ", &
    report%minimum_gap_norm, " at ", report%minimum_gap_x, &
    report%minimum_gap_y

  if (argument_count == 2) then
    call get_command_argument(2, output_path)
    call write_spinful_field_map_2d( &
      trim(output_path), mesh, state, "legacy_split_2d_import", &
      "legacy_supplied_field")
    print '(a,a)', "wrote FermiForge field map: ", trim(output_path)
  end if
end program import_legacy_split_field_map_2d
