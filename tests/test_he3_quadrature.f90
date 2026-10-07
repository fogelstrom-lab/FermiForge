program test_he3_quadrature
  use he3_kinds, only : rk
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t, &
                             read_legacy_gauss_table, &
                             read_legacy_ozaki_table, generate_ozaki_quadrature
  implicit none

  type(angular_quadrature_3d_t) :: angular
  type(ozaki_quadrature_t) :: energy
  type(ozaki_quadrature_t) :: generated
  character(len=1024) :: gauss_path, ozaki_path
  real(rk) :: expected_first(3), moment(3), phi, planar_magnitude

  if (command_argument_count() == 1) then
    call get_command_argument(1, ozaki_path)
    select case(trim(ozaki_path))
    case('zero_temperature')
      call generate_ozaki_quadrature(0.0_rk,50.0_rk,generated)
    case('negative_cutoff')
      call generate_ozaki_quadrature(0.3_rk,-1.0_rk,generated)
    case('excess_cutoff')
      call generate_ozaki_quadrature(0.3_rk,1.e6_rk,generated)
    end select
    stop 0
  end if

  if (command_argument_count() /= 2) &
    error stop "test_he3_quadrature needs Gauss and Ozaki table paths"
  call get_command_argument(1, gauss_path)
  call get_command_argument(2, ozaki_path)

  call read_legacy_gauss_table(trim(gauss_path), 48, 11, angular)
  call require(angular%is_valid(), "legacy angular quadrature is invalid")
  call require(angular%direction_count() == 528, &
               "legacy angular direction count changed")
  call require(abs(sum(angular%weight) - 1.0_rk) < 2.0e-13_rk, &
               "legacy angular weights are not normalized")
  moment = angular%second_moment()
  call require(maxval(abs(moment - 1.0_rk / 3.0_rk)) < 3.0e-13_rk, &
               "legacy angular quadrature has the wrong second moment")

  phi = 0.5_rk * acos(-1.0_rk) / 48.0_rk
  planar_magnitude = sqrt(1.0_rk - (-0.9782286581461_rk)**2)
  expected_first = [cos(phi) * planar_magnitude, &
                    sin(phi) * planar_magnitude, -0.9782286581461_rk]
  call require(maxval(abs(angular%momentum(:, 1) - expected_first)) < &
               3.0e-15_rk, "legacy angular direction ordering changed")
  call require(abs(angular%weight(1) - &
               0.5_rk * 0.0556685671161737_rk / 48.0_rk) < 2.0e-18_rk, &
               "legacy first angular weight changed")

  call read_legacy_ozaki_table(trim(ozaki_path), energy)
  call generate_ozaki_quadrature(0.3_rk,50.0_rk,generated)
  call require(generated%pole_count()==energy%pole_count(), "generated count differs")
  call require(maxval(abs(generated%pole-energy%pole))<1.e-10_rk, "generated poles differ")
  call require(maxval(abs(generated%residue-energy%residue))<1.e-9_rk, "generated residues differ")
  call require(abs(generated%gap_prefactor-energy%gap_prefactor)<1.e-12_rk, "generated gap prefactor differs")
  call generate_ozaki_quadrature(0.1_rk,50.0_rk,generated)
  call require(generated%is_valid().and.generated%pole_count()==14, "low temperature generation")
  call require(energy%is_valid(), "legacy Ozaki quadrature is invalid")
  call require(energy%pole_count() == 8, "legacy Ozaki pole count changed")
  call require(abs(energy%temperature - 0.3_rk) < tiny(1.0_rk), &
               "legacy Ozaki temperature changed")
  call require(energy%cutoff_index == 50, "legacy Ozaki cutoff changed")
  call require(abs(energy%pole(1) - 0.15_rk) < 2.0e-16_rk .and. &
               abs(energy%residue(8) - 54.95372645203803_rk) < 2.0e-13_rk, &
               "legacy Ozaki table values changed")
  call require(energy%gap_prefactor > 0.0_rk .and. &
               energy%gap_prefactor < 1.0_rk, &
               "legacy Ozaki gap prefactor is outside its expected range")

  print '(a)', "3He quadrature tests passed"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_he3_quadrature
