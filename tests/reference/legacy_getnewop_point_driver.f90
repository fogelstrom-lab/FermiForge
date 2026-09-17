program legacy_getnewop_point_driver
  use Global_variables
  use Initialisation, only : init_calc
  use Interpolations, only : intord_v
  use riccati, only : makeprops
  implicit none

  complex :: energy, self_energy(4, -mx:mx), mapped_value(12)
  character(len=32) :: argument
  integer :: component, direction, polar_node, pole, radial_index

  radial_index = 24
  if (command_argument_count() >= 1) then
    call get_command_argument(1, argument)
    read(argument, *) radial_index
  end if
  if (radial_index < 0 .or. radial_index > nx) &
    error stop "radial reference point lies outside the source grid"

  ! init_calc deliberately reads the unmodified new_src input and data files.
  call init_calc
  mapped_value = cmplx(0.0, 0.0)

  do direction = 1, tmax
    do polar_node = 1, 11
      call intord_v(direction, polar_node, radial_index, self_energy)
      do pole = 1, Ncmax
        energy = cmplx(zp(pole), 0.0)
        call makeprops(direction, polar_node, pole, energy, self_energy, &
                       mapped_value)
      end do
    end do
  end do
  mapped_value(1:9) = esum * mapped_value(1:9)

  write(*, '(a,i0)') "REFERENCE_RADIAL_INDEX ", radial_index
  write(*, '(a,es25.17)') "REFERENCE_RADIUS ", xgrid(radial_index)
  do component = 1, 12
    write(*, '(a,1x,i2,2(1x,es25.17))') "REFERENCE_COMPONENT", &
      component, real(mapped_value(component)), aimag(mapped_value(component))
  end do
end program legacy_getnewop_point_driver
