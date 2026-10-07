module radial_iteration_diagnostics
  use he3_kinds, only: rk
  use spinful_state_2d, only: spinful_state_2d_t
  use new_src_iteration_layout_2d, only: pack_new_src_masked_iteration_vector_2d
  use barzilai_borwein_mixing, only: bb_mixing_t
  implicit none
  private
  public :: write_radial_iteration_diagnostics
contains
  subroutine write_radial_iteration_diagnostics(prefix,k,state,mapped,next,active,bb,alpha,inner,outer,engine)
    character(len=*),intent(in) :: prefix
    character(len=*),intent(in) :: engine
    integer,intent(in) :: k
    type(spinful_state_2d_t),intent(in) :: state,mapped,next
    logical,intent(in) :: active(:)
    type(bb_mixing_t),intent(in) :: bb
    real(rk),intent(in) :: alpha,inner,outer
    real(rk),allocatable :: x(:),f(:),xn(:),r(:)
    real(rk) :: squares(3),maxima(3),radius, value
    integer :: counts(3),u,i,p,j,n,zone
    character(len=6) :: suffix
    if(.not.allocated(state%radial_point)) error stop 'signed diagnostics require radial symmetry'
    n=count(active)
    if(n/=size(state%radial_point)) error stop 'diagnostic radial mask mismatch'
    call pack_new_src_masked_iteration_vector_2d(state,active,x)
    call pack_new_src_masked_iteration_vector_2d(mapped,active,f)
    call pack_new_src_masked_iteration_vector_2d(next,active,xn)
    r=f-x
    write(suffix,'(i6.6)') k
    open(newunit=u,file=prefix//'.signed.map'//suffix//'.dat',status='replace',action='write')
    write(u,'(a)') '# packed_index state residual_F_minus_X applied_update_next_minus_X'
    do i=1,size(x)
      write(u,'(i8,3es25.16e3)') i,x(i),r(i),xn(i)-x(i)
    end do
    close(u)
    if(k==1) then
      open(newunit=u,file=prefix//'.signed.layout.dat',status='replace',action='write')
      write(u,'(a)') '# block-major: Re Axx, Im Axx, Re Axy, Im Axy, ... Ayx ... Azz, mf_x,mf_y,mf_z'
      write(u,'(a)') '# each of 21 blocks has N rows; index=(component-1)*N+radial_index'
      write(u,'(a)') '# radial_index mesh_point radius zone; unweighted zones r<=inner, inner<r<=outer, r>outer'
      do p=1,n
        radius=state%radial_coordinate(p)
        zone=1
        if(radius>inner) zone=2
        if(radius>outer) zone=3
        write(u,'(2i8,es25.16e3,i4)') p,state%radial_point(p),radius,zone
      end do
      close(u)
    end if
    squares=0; maxima=0; counts=0
    do p=1,n
      radius=state%radial_coordinate(p)
      zone=1
      if(radius>inner) zone=2
      if(radius>outer) zone=3
      do j=1,21
        value=r((j-1)*n+p)
        squares(zone)=squares(zone)+value**2
        maxima(zone)=max(maxima(zone),abs(value))
        counts(zone)=counts(zone)+1
      end do
    end do
    if(k==1) then
      open(newunit=u,file=prefix//'.regional_diagnostics.dat',status='replace',action='write')
      write(u,'(a)') '# iteration step_parameter core_rms transition_rms outer_rms core_max transition_max outer_max'
    else
      open(newunit=u,file=prefix//'.regional_diagnostics.dat',status='old',position='append',action='write')
    end if
    write(u,'(i8,7es25.16e3)') k,alpha,sqrt(squares/real(max(counts,1),rk)),maxima
    close(u)
    ! BB secant metadata has no meaning for Anderson or Polyak. Common signed
    ! vectors allow offline reconstruction of observed secants for any engine.
    if(engine/='bb') return
    if(k==1) then
      open(newunit=u,file=prefix//'.bb_diagnostics.dat',status='replace',action='write')
      write(u,'(a)') '# iteration alpha formula status limited secant_valid alignment_valid ' // &
        'signed_cosine residual_alignment s_norm y_norm raw_BB1 raw_BB2 ' // &
        'core_rms transition_rms outer_rms core_max transition_max outer_max'
    else
      open(newunit=u,file=prefix//'.bb_diagnostics.dat',status='old',position='append',action='write')
    end if
    write(u,'(i8,es25.16e3,5i4,12es25.16e3)') k,alpha,bb%formula,bb%status, &
      merge(1,0,bb%limited),merge(1,0,bb%secant_valid),merge(1,0,bb%alignment_valid), &
      bb%signed_cosine,bb%residual_alignment,bb%step_norm,bb%secant_norm,bb%raw_bb1,bb%raw_bb2, &
      sqrt(squares/real(max(counts,1),rk)),maxima
    close(u)
  end subroutine
end module
