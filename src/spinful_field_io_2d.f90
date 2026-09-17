module spinful_field_io_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t
  implicit none
  private

  public :: write_spinful_field_map_2d

contains

  subroutine write_spinful_field_map_2d( &
      path, mesh, state, source_label, endpoint_policy)
    character(len=*), intent(in) :: path
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    character(len=*), intent(in) :: source_label, endpoint_policy

    character(len=512) :: output_message
    integer :: output_status, output_unit

    if (.not. state%is_valid_for(mesh)) &
      error stop "cannot write a 2D state that does not match its mesh"
    open(newunit=output_unit, file=trim(path), status="replace", &
         action="write", iostat=output_status, iomsg=output_message)
    if (output_status /= 0) &
      error stop "cannot open 2D field output: " // trim(output_message)
    call write_field_map_records( &
      output_unit, mesh, state, source_label, endpoint_policy)
    close(output_unit)
  end subroutine write_spinful_field_map_2d


  subroutine write_field_map_records( &
      unit, mesh, state, source_label, endpoint_policy)
    integer, intent(in) :: unit
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    character(len=*), intent(in) :: source_label, endpoint_policy

    real(rk) :: pair_density, row(24), x, y
    integer :: orbital, point, position, spin, x_node, y_node

    write(unit, '(a)') "# FermiForge 2D field map version=1"
    write(unit, '(a,a)') "# source=", trim(source_label)
    write(unit, '(a)') "# density_kind=pair_density"
    write(unit, '(a)') "# current_kind=current_related_mean_field"
    write(unit, '(a,a)') "# endpoint_policy=", trim(endpoint_policy)
    write(unit, '(a,a)') "# mesh_kind=", &
      merge("uniform    ", "rectilinear", mesh%is_uniform())
    write(unit, '(a,es24.16e3)') "# minimum_spacing=", &
      mesh%minimum_spacing()
    write(unit, '(a,es24.16e3)') "# maximum_spacing=", &
      mesh%maximum_spacing()
    write(unit, '(a,i0)') "# nx_points=", mesh%x_point_count()
    write(unit, '(a,i0)') "# ny_points=", mesh%y_point_count()
    write(unit, '(a)') &
      "# columns=x y A_xx_re A_xx_im A_xy_re A_xy_im A_xz_re A_xz_im " // &
      "A_yx_re A_yx_im A_yy_re A_yy_im A_yz_re A_yz_im " // &
      "A_zx_re A_zx_im A_zy_re A_zy_im A_zz_re A_zz_im " // &
      "pair_density j_x j_y j_z"
    do y_node = 0, mesh%y_cell_count()
      y = mesh%y_coordinate(y_node)
      do x_node = 0, mesh%x_cell_count()
        x = mesh%x_coordinate(x_node)
        point = mesh%point_index(x_node, y_node)
        row = 0.0_rk
        row(1:2) = [x, y]
        position = 2
        do spin = 1, 3
          do orbital = 1, 3
            position = position + 1
            row(position) = real(state%order_parameter(spin, orbital, point), rk)
            position = position + 1
            row(position) = aimag(state%order_parameter(spin, orbital, point))
          end do
        end do
        pair_density = sum(abs(state%order_parameter(:, :, point))**2) / &
          3.0_rk
        row(21) = pair_density
        row(22:24) = state%current_mean_field(:, point)
        write(unit, '(*(es24.16e3,1x))') row
      end do
    end do
  end subroutine write_field_map_records

end module spinful_field_io_2d
