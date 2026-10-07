program test_bulk_gap
  use he3_kinds, only: rk
  use he3_fermi_liquid, only: resolve_fs1
  use he3_quadrature, only: ozaki_quadrature_t, generate_ozaki_quadrature
  use he3_bulk_gap, only: bulk_gap_legacy,bulk_gap_ozaki,resolve_bulk_gap
  implicit none
  type(ozaki_quadrature_t) :: q
  real(rk) :: gap, previous, t, residual
  integer :: i
  character(len=16) :: mode
  previous=1.0_rk
  gap=99.0_rk
  call resolve_fs1(5.4_rk,gap)
  if(abs(gap-1.9285714285714286_rk)>1.e-14_rk) error stop 'Fs1 conversion'
  call resolve_fs1(huge(1.0_rk),gap)
  if(abs(gap-1.9285714285714286_rk)>1.e-14_rk) error stop 'legacy feedback preservation'
  call resolve_fs1(0.0_rk,gap)
  if(gap/=0.0_rk) error stop 'zero Fs1'
  ! Original new_src/gap.f90 prints 0.279946 at T=.30, irep=1.
  if(abs(bulk_gap_legacy(0.3_rk)-0.279946_rk)>0.5e-6_rk) error stop 'original bulk-gap reference'
  do i=1,9
    t=real(i,rk)/10.0_rk
    call generate_ozaki_quadrature(t,50.0_rk,q)
    gap=bulk_gap_ozaki(q)
    if(gap<=0.or.gap>=previous) error stop 'bulk gap temperature monotonicity'
    residual=gap-q%gap_prefactor*sum(q%residue*gap/sqrt(q%pole**2+gap**2))
    if(abs(residual)>1.e-13_rk) error stop 'bulk gap does not satisfy discrete solver equation'
    if(abs(gap-bulk_gap_legacy(t))>1.e-4_rk) error stop 'legacy gap differs unexpectedly'
    previous=gap
  end do
  if(bulk_gap_legacy(1.0_rk)/=0.0_rk) error stop 'gap at Tc'
  call generate_ozaki_quadrature(1.0_rk,50.0_rk,q)
  if(bulk_gap_ozaki(q)/=0.0_rk) error stop 'Ozaki gap at Tc'
  call generate_ozaki_quadrature(0.3_rk,50.0_rk,q)
  mode='auto'
  gap=0.123_rk
  call resolve_bulk_gap(mode,'file',q,gap)
  if(mode/='manual'.or.gap/=0.123_rk) error stop 'legacy inputs changed'
  mode='auto'
  call resolve_bulk_gap(mode,'generate',q,gap)
  if(mode/='ozaki'.or.abs(gap-bulk_gap_ozaki(q))>1.e-14_rk) error stop 'automatic gap selection'
  print *, 'T=0.3 bulk gaps [legacy, Ozaki]: ',bulk_gap_legacy(0.3_rk),gap
end program
