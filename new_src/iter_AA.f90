module anderson_iteration_adapter
  use, intrinsic :: iso_fortran_env, only : real64
  use global_variables, only : aa_pmax, iiter, nop, vecl
  use mpi_variables, only : myid
  use legacy_anderson_mixing, only : legacy_anderson_t, anderson_report_t
  use newses, only : getnewop
  use packses, only : pacvector
  implicit none
  private

  integer, parameter :: history_limit = 10
  real(real64), parameter :: progress_threshold = 0.1_real64

  public :: iterateAA

contains

  subroutine iterateAA(tolA, tolL, p, md)
    ! Compatibility adapter between the experimental free-form vortex driver
    ! and the modern translation of the original citerat algorithm in iter.f.
    ! The nonlinear map is still evaluated collectively by getnewop; only rank
    ! zero owns and updates the Anderson history.
    real, intent(out) :: tolA, tolL, p
    integer, intent(out) :: md

    type(legacy_anderson_t), save :: accelerator
    type(anderson_report_t) :: report
    integer, save :: lower(nop), upper(nop)
    real(real64), allocatable, save :: current(:), mapped(:), next(:), residual(:)
    logical, save :: initialized = .false.
    integer :: i, number_of_values
    real(real64) :: current_norm, residual_norm

    ! Values on non-root ranks are replaced by broadcasts in qcv immediately
    ! after this routine returns. Initializing them here also makes every
    ! intent(out) argument well-defined.
    tolA = huge(1.0)
    tolL = huge(1.0)
    p = 0.0
    md = 0

    if (myid == 0) then
      if (iiter == 1 .or. .not. initialized) then
        lower(1) = 1
        upper(1) = vecl
        do i = 2, nop
          lower(i) = lower(i - 1) + vecl
          upper(i) = upper(i - 1) + vecl
        end do
        number_of_values = upper(nop)

        if (allocated(current)) deallocate(current, mapped, next, residual)
        allocate(current(number_of_values), mapped(number_of_values), &
                 next(number_of_values), residual(number_of_values))

        call accelerator%initialize(number_of_values, history_limit, &
                                    progress_threshold, real(aa_pmax, real64))
        initialized = .true.
      end if

      call pacvector(lower, upper, current, 1)
    end if

    ! getnewop evaluates G(current) on all MPI ranks and reduces the mapped
    ! fields to rank zero.
    call getnewop

    if (myid == 0) then
      ! pacvector mode 2 replaces its input by G(current)-current. Supplying a
      ! copy of current therefore recovers both the residual and mapped state.
      residual = current
      call pacvector(lower, upper, residual, 2)
      mapped = current + residual

      ! qcv owns the convergence test.  A zero tolerance keeps the accelerator
      ! from applying its legacy absolute-error stopping test while qcv uses
      ! the relative diagnostics below, as do the NN, BR, and BB adapters.
      call accelerator%update(current, mapped, 0.0_real64, next, report)
      call pacvector(lower, upper, next, 3)

      ! Retain the two convergence diagnostics used by the other new_src
      ! iteration schemes: a global relative 2-norm and a scaled largest-
      ! component residual.  They were inadvertently made identical in the
      ! first adapter revision.
      current_norm = sqrt(max(dot_product(current, current), 0.0_real64))
      if (current_norm <= tiny(current_norm)) current_norm = 1.0_real64
      residual_norm = sqrt(max(dot_product(residual, residual), 0.0_real64))
      tolA = real(residual_norm / current_norm, kind(tolA))
      tolL = real(real(size(residual), real64) * &
                  report%max_residual / current_norm, kind(tolL))
      p = real(report%mixing, kind(p))
      md = report%history_size
    end if
  end subroutine iterateAA

end module anderson_iteration_adapter
