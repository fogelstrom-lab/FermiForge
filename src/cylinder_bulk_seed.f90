module cylinder_bulk_seed
  use he3_kinds, only: rk
  implicit none
  private
  public :: make_cylinder_bulk_seed
  public :: make_cylinder_texture_seed
contains
  subroutine make_cylinder_texture_seed(texture,x,y,radius,amplitude,seed)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(len=*), intent(in) :: texture
    real(rk), intent(in) :: x,y,radius,amplitude
    complex(rk), intent(out) :: seed(3,3)
    real(rk) :: r,phi,beta,c,s,eta
    complex(rk) :: phase,orbital(3)
    if(.not.all(ieee_is_finite([x,y,radius,amplitude]))) error stop 'nonfinite texture input'
    if(radius<=0.or.amplitude<=0) error stop 'invalid texture scale'
    r=hypot(x,y)
    if(r>radius+128*epsilon(radius)*max(1._rk,radius)) error stop 'texture point outside cylinder'
    seed=0
    if(texture/='a_mermin_ho'.and.texture/='a_planar'.and.texture/='a_panam') &
      error stop 'unknown cylinder texture'
    if(r==0) then
      if(texture=='a_panam') then
        seed(1,1)=amplitude; seed(1,2)=cmplx(0._rk,amplitude,rk)
      else
        seed(3,1)=amplitude; seed(3,2)=cmplx(0._rk,amplitude,rk)
      end if
      return
    end if
    phi=atan2(y,x); beta=acos(-1._rk)*min(r/radius,1._rk)/2
    c=cos(beta); s=sin(beta)
    if(texture=='a_mermin_ho'.or.texture=='a_panam') then
      ! Smooth pure A: d=z, l=sin(beta)e_r+cos(beta)e_z.
      ! e^(i phi) cancels the triad's coordinate singularity at the origin.
      ! One circulation quantum at the wall; no extra vortex multiplier.
      phase=cmplx(x/r,y/r,rk)
      seed(3,1)=amplitude*phase*cmplx(c*cos(phi),-sin(phi),rk)
      seed(3,2)=amplitude*phase*cmplx(c*sin(phi),cos(phi),rk)
      seed(3,3)=-amplitude*phase*s
      if(texture=='a_panam') then
        ! Smooth schematic hyperbolic-like spin director, NOT a digitization
        ! or GL minimizer of Takagi Fig.2. The orbital l texture stays MH.
        eta=-acos(-1._rk)*x*y/(2*radius**2)
        orbital=seed(3,:); seed=0
        seed(1,:)=cos(eta)*orbital; seed(2,:)=sin(eta)*orbital
      end if
    else
      ! Literal formula from FermiForge-notes.pdf, with one global scale.
      ! NOT the conventional planar phase: mixed A/polar trial state, with
      ! a polar point at (x,y)=(0,R/2). No local normalization/correction.
      seed(3,1)=amplitude*(c-sin(phi)*s)
      seed(3,2)=amplitude*cmplx(cos(phi)*s,c,rk)
      seed(3,3)=cmplx(0._rk,amplitude*s,rk)
    end if
  end subroutine

  subroutine make_cylinder_bulk_seed(initialization,physical_winding,gap,seed,symmetry_winding)
    integer, intent(in) :: initialization
    real(rk), intent(in) :: physical_winding,gap
    complex(rk), intent(out) :: seed(3,3)
    real(rk), intent(out) :: symmetry_winding
    integer :: j
    if(physical_winding/=0) error stop 'bulk cylinder seeds require zero physical phase winding'
    seed=cmplx(0._rk,0._rk,rk)
    select case(initialization)
    case(-1)
      do j=1,3
        seed(j,j)=gap
      end do
      symmetry_winding=0
    case(-2)
      ! new_src bulkA: d_z (p_x+i p_y), with its original seed amplitude.
      seed(3,1)=sqrt(2._rk)*gap
      seed(3,2)=cmplx(0._rk,sqrt(2._rk)*gap,rk)
      ! Only C_(0,+) is nonzero. In exp[i(m-s-k)phi], m=1
      ! makes its Cartesian seed constant: intrinsic chirality, NOT a vortex.
      symmetry_winding=1
    case default
      error stop 'cylinder initialization must be -1 (B) or -2 (A)'
    end select
  end subroutine
end module
