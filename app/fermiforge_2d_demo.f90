program fermiforge_2d_demo
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use bilinear_field_sampler_2d, only : sample_spinful_state_2d
  implicit none

  character(len=512) :: input_file, output_file, output_override
  character(len=512) :: input_message
  integer :: input_status, input_unit, output_unit
  integer :: number_of_x_cells, number_of_y_cells, phase_winding
  real(rk) :: x_minimum, x_maximum, y_minimum, y_maximum
  real(rk) :: core_width, texture_width, off_diagonal_scale
  real(rk) :: current_scale, axial_current_scale
  logical :: output_at_cell_centers

  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state

  namelist /demo_2d/ x_minimum, x_maximum, number_of_x_cells, &
                     y_minimum, y_maximum, number_of_y_cells, &
                     core_width, texture_width, off_diagonal_scale, &
                     current_scale, axial_current_scale, phase_winding, &
                     output_at_cell_centers, output_file

  x_minimum = -8.0_rk
  x_maximum = 8.0_rk
  number_of_x_cells = 128
  y_minimum = -8.0_rk
  y_maximum = 8.0_rk
  number_of_y_cells = 128
  core_width = 1.0_rk
  texture_width = 3.0_rk
  off_diagonal_scale = 0.18_rk
  current_scale = 0.40_rk
  axial_current_scale = 0.035_rk
  phase_winding = 1
  output_at_cell_centers = .true.
  output_file = "fermiforge_2d_demo.dat"

  call get_command_argument(1, input_file)
  if (len_trim(input_file) == 0) input_file = "examples/2d_sampler_demo.nml"

  open(newunit=input_unit, file=trim(input_file), status="old", action="read", &
       iostat=input_status, iomsg=input_message)
  if (input_status /= 0) &
    error stop "cannot open 2D demo input: " // trim(input_message)
  read(input_unit, nml=demo_2d, iostat=input_status, iomsg=input_message)
  close(input_unit)
  if (input_status /= 0) &
    error stop "cannot read 2D demo input: " // trim(input_message)

  call get_command_argument(2, output_override)
  if (len_trim(output_override) > 0) output_file = output_override

  call validate_input()
  call make_uniform_cartesian_mesh(x_minimum, x_maximum, number_of_x_cells, &
                                   y_minimum, y_maximum, number_of_y_cells, mesh)
  call allocate_spinful_state_2d(mesh, state)
  call initialize_manufactured_vortex(mesh, state)

  open(newunit=output_unit, file=trim(output_file), status="replace", &
       action="write", iostat=input_status, iomsg=input_message)
  if (input_status /= 0) &
    error stop "cannot open 2D demo output: " // trim(input_message)
  call write_field_map(output_unit, mesh, state)
  close(output_unit)

  write(*, '(a)') "FermiForge manufactured 2D field map"
  write(*, '(a,a)') "input:  ", trim(input_file)
  write(*, '(a,a)') "output: ", trim(output_file)

contains

  subroutine validate_input()
    if (x_maximum <= x_minimum .or. y_maximum <= y_minimum) &
      error stop "2D demo coordinate bounds must be strictly increasing"
    if (number_of_x_cells < 2 .or. number_of_y_cells < 2) &
      error stop "2D demo needs at least two cells in each direction"
    if (core_width <= 0.0_rk .or. texture_width <= 0.0_rk) &
      error stop "2D demo core and texture widths must be positive"
    if (off_diagonal_scale < 0.0_rk .or. current_scale < 0.0_rk .or. &
        axial_current_scale < 0.0_rk) &
      error stop "2D demo amplitudes cannot be negative"
    if (len_trim(output_file) == 0) error stop "2D demo output filename is empty"
  end subroutine validate_input


  subroutine initialize_manufactured_vortex(field_mesh, field_state)
    type(cartesian_mesh_2d_t), intent(in) :: field_mesh
    type(spinful_state_2d_t), intent(inout) :: field_state

    complex(rk) :: angular_phase, gauge_phase, vortex_amplitude
    real(rk) :: angle, core_envelope, radius, radius_squared
    integer :: orbital, point, spin, x_node, y_node
    real(rk) :: x, y

    do y_node = 0, field_mesh%y_cell_count()
      y = field_mesh%y_coordinate(y_node)
      do x_node = 0, field_mesh%x_cell_count()
        x = field_mesh%x_coordinate(x_node)
        point = field_mesh%point_index(x_node, y_node)
        radius_squared = x * x + y * y
        radius = sqrt(radius_squared)
        if (radius > tiny(radius)) then
          angle = atan2(y, x)
        else
          angle = 0.0_rk
        end if

        gauge_phase = cmplx(cos(real(phase_winding, rk) * angle), &
                            sin(real(phase_winding, rk) * angle), kind=rk)
        vortex_amplitude = cmplx(tanh(radius / core_width), 0.0_rk, kind=rk) * &
                           gauge_phase
        core_envelope = exp(-radius_squared / (texture_width * texture_width))

        do spin = 1, 3
          do orbital = 1, 3
            if (spin == orbital) then
              field_state%order_parameter(spin, orbital, point) = &
                vortex_amplitude * cmplx( &
                1.0_rk + 0.06_rk * real(spin - 2, rk) * core_envelope, &
                0.0_rk, kind=rk)
            else
              angular_phase = cmplx( &
                cos(real(spin - orbital, rk) * angle), &
                sin(real(spin - orbital, rk) * angle), kind=rk)
              field_state%order_parameter(spin, orbital, point) = &
                cmplx(off_diagonal_scale * core_envelope, 0.0_rk, kind=rk) * &
                vortex_amplitude * &
                angular_phase
            end if
          end do
        end do

        field_state%current_mean_field(1, point) = &
          -current_scale * y / (radius_squared + core_width * core_width)
        field_state%current_mean_field(2, point) = &
           current_scale * x / (radius_squared + core_width * core_width)
        field_state%current_mean_field(3, point) = &
          axial_current_scale * exp(-radius_squared / (core_width * core_width))
      end do
    end do
  end subroutine initialize_manufactured_vortex


  subroutine write_field_map(unit, field_mesh, field_state)
    integer, intent(in) :: unit
    type(cartesian_mesh_2d_t), intent(in) :: field_mesh
    type(spinful_state_2d_t), intent(in) :: field_state

    complex(rk) :: sampled_order(3, 3)
    real(rk) :: sampled_current(3), row(24)
    real(rk) :: pair_density, x, y
    integer :: number_of_output_x_points, number_of_output_y_points
    integer :: orbital, output_x, output_y, position, spin
    logical :: inside

    if (output_at_cell_centers) then
      number_of_output_x_points = field_mesh%x_cell_count()
      number_of_output_y_points = field_mesh%y_cell_count()
    else
      number_of_output_x_points = field_mesh%x_point_count()
      number_of_output_y_points = field_mesh%y_point_count()
    end if

    write(unit, '(a)') "# FermiForge 2D field map version=1"
    write(unit, '(a)') "# source=manufactured_axial_vortex"
    write(unit, '(a)') "# density_kind=pair_density"
    write(unit, '(a)') "# current_kind=manufactured_current"
    write(unit, '(a,i0)') "# phase_winding=", phase_winding
    write(unit, '(a,i0)') "# nx_points=", number_of_output_x_points
    write(unit, '(a,i0)') "# ny_points=", number_of_output_y_points
    write(unit, '(a)') &
      "# columns=x y A_xx_re A_xx_im A_xy_re A_xy_im A_xz_re A_xz_im " // &
      "A_yx_re A_yx_im A_yy_re A_yy_im A_yz_re A_yz_im " // &
      "A_zx_re A_zx_im A_zy_re A_zy_im A_zz_re A_zz_im " // &
      "pair_density j_x j_y j_z"

    do output_y = 0, number_of_output_y_points - 1
      if (output_at_cell_centers) then
        y = field_mesh%y_coordinate(output_y) + 0.5_rk * field_mesh%y_spacing()
      else
        y = field_mesh%y_coordinate(output_y)
      end if
      do output_x = 0, number_of_output_x_points - 1
        if (output_at_cell_centers) then
          x = field_mesh%x_coordinate(output_x) + 0.5_rk * field_mesh%x_spacing()
        else
          x = field_mesh%x_coordinate(output_x)
        end if

        call sample_spinful_state_2d(field_mesh, field_state, x, y, &
                                     sampled_order, sampled_current, inside)
        if (.not. inside) error stop "2D demo attempted an outside field sample"

        row = 0.0_rk
        row(1:2) = [x, y]
        position = 2
        do spin = 1, 3
          do orbital = 1, 3
            position = position + 1
            row(position) = real(sampled_order(spin, orbital), kind=rk)
            position = position + 1
            row(position) = aimag(sampled_order(spin, orbital))
          end do
        end do
        pair_density = sum(abs(sampled_order)**2) / 3.0_rk
        row(21) = pair_density
        row(22:24) = sampled_current
        write(unit, '(*(es24.16e3,1x))') row
      end do
    end do
  end subroutine write_field_map

end program fermiforge_2d_demo
