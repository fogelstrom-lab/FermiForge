! Isotropic (irep=1, zero-flow) reduction of new_src/bulkgap.f90.
! Energies and gap in the existing solver units, Delta/(2*pi*kB*Tc).
module he3_bulk_gap
  use he3_kinds, only: rk
  use he3_quadrature, only: ozaki_quadrature_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: bulk_gap_ozaki, bulk_gap_legacy, resolve_bulk_gap
contains
  function solve_gap(t,pole,weight) result(gap)
    real(rk), intent(in) :: t,pole(:),weight(:)
    real(rk) :: gap,lo,hi,mid,value
    integer :: iteration
    if (.not.ieee_is_finite(t).or.t<=0.0_rk.or.t>1.0_rk) &
      error stop 'bulk gap requires 0 < T/Tc <= 1'
    if(t==1.0_rk) then
      gap=0.0_rk
      return
    end if
    lo=0.0_rk
    hi=1.0_rk
    do iteration=1,60
      if(residual(hi)>0.0_rk) exit
      hi=2.0_rk*hi
    end do
    if(iteration>60) error stop 'cannot bracket bulk gap'
    do iteration=1,200
      mid=0.5_rk*(lo+hi)
      value=residual(mid)
      if(value>0.0_rk) then
        hi=mid
      else
        lo=mid
      end if
      if(hi-lo<2.e-14_rk) exit
    end do
    gap=0.5_rk*(lo+hi)
  contains
    function residual(d) result(f)
      real(rk), intent(in) :: d
      real(rk) :: f,e(size(pole))
      e=sqrt(pole**2+d**2)
      ! Stable form of log(T) + T sum w*(1/pole - 1/E).
      f=log(t)+t*sum(weight*d**2/(pole*e*(e+pole)))
    end function
  end function

  function bulk_gap_ozaki(energy) result(gap)
    type(ozaki_quadrature_t), intent(in) :: energy
    real(rk) :: gap
    if(.not.energy%is_valid()) error stop 'invalid bulk-gap quadrature'
    gap=solve_gap(energy%temperature,energy%pole,energy%residue)
  end function

  function bulk_gap_legacy(t) result(gap)
    real(rk), intent(in) :: t
    real(rk) :: gap
    real(rk), allocatable :: poles(:),weights(:)
    integer :: n,i
    if(.not.ieee_is_finite(t).or.t<=0.0_rk.or.t>1.0_rk) &
      error stop 'legacy bulk gap requires 0 < T/Tc <= 1'
    if(t<15.0_rk/1000000.0_rk) error stop 'legacy bulk-gap sum exceeds supported size'
    n=int(15.0_rk/t+0.00001_rk)+1
    allocate(poles(n),weights(n))
    do i=1,n
      poles(i)=t*(real(i-1,rk)+0.5_rk)
    end do
    weights=1.0_rk
    gap=solve_gap(t,poles,weights)
  end function

  subroutine resolve_bulk_gap(mode,ozaki_mode,energy,gap)
    character(len=*), intent(inout) :: mode
    character(len=*), intent(in) :: ozaki_mode
    type(ozaki_quadrature_t), intent(in) :: energy
    real(rk), intent(inout) :: gap
    if(mode=='auto') then
      if(ozaki_mode=='generate') then
        mode='ozaki'
      else
        mode='manual'
      end if
    end if
    select case(trim(mode))
    case('ozaki')
      gap=bulk_gap_ozaki(energy)
    case('legacy')
      gap=bulk_gap_legacy(energy%temperature)
    case('manual')
    case default
      error stop 'bulk_gap_mode must be auto, ozaki, legacy or manual'
    end select
    if(.not.ieee_is_finite(gap).or.gap<=0.0_rk) &
      error stop 'vortex bulk gap must be finite and positive (T < Tc)'
  end subroutine
end module
