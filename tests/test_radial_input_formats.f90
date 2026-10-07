program test_radial_input_formats
  use Initialisation, only : read_radial_input
  use Global_variables, only : vort,icyl,aa0,tmax,istart,ittyp,itmax,aa_pmax,Rx
  implicit none
  real :: temperature,tolerance,extent
  call read_radial_input(temperature,tolerance,extent)
  write(*,'(12(es24.16,1x))') temperature,vort,real(icyl),extent,Rx,aa0, &
    real(tmax),tolerance,real(istart),real(ittyp),real(itmax),aa_pmax
end program
