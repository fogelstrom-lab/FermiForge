module mpi_shutdown_audit
  use mpi_f08
  use iso_fortran_env, only: int64, real64
  implicit none
contains
  subroutine finalize_with_audit(prefix, ierr)
    character(len=*), intent(in) :: prefix
    integer, intent(out) :: ierr
    integer :: rank, code, unit, ios, length
    logical :: finalized
    character(len=32) :: suffix
    character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: version
    integer(int64) :: start_tick,end_tick,tick_rate
    character(len=128) :: setting
    call MPI_Comm_rank(MPI_COMM_WORLD,rank,code)
    if(code /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,code,ierr)
    call MPI_Get_library_version(version,length,code)
    write(suffix,'(".shutdown.rank",i6.6,".txt")') rank
    open(newunit=unit,file=trim(prefix)//trim(suffix),status='replace',iostat=ios)
    if(ios /= 0) then
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      error stop 'cannot open MPI shutdown audit'
    end if
    if(code == MPI_SUCCESS) write(unit,'(a)') trim(version(:length))
    call get_environment_variable('PMIX_MCA_pmix_finalize_timeout',setting,status=code)
    if(code==0) write(unit,'(a,a)') 'PMIX_MCA_pmix_finalize_timeout=',trim(setting)
    call get_environment_variable('PMIX_MCA_pmix_client_base_verbose',setting,status=code)
    if(code==0) write(unit,'(a,a)') 'PMIX_MCA_pmix_client_base_verbose=',trim(setting)
    write(unit,'(a)') 'entering_MPI_Finalize'
    flush(unit)
    call system_clock(start_tick,tick_rate)
    call MPI_Finalize(ierr)
    call system_clock(end_tick)
    write(unit,'(a,i0)') 'MPI_Finalize_return_code=',ierr
    if(tick_rate>0) write(unit,'(a,f14.6)') 'MPI_Finalize_elapsed_seconds=', &
      real(end_tick-start_tick,real64)/real(tick_rate,real64)
    ! MPI_Finalized is explicitly valid after finalization.
    call MPI_Finalized(finalized,code)
    write(unit,'(a,l1)') 'MPI_Finalized=',finalized
    write(unit,'(a,i0)') 'MPI_Finalized_return_code=',code
    close(unit)
  end subroutine finalize_with_audit
end module mpi_shutdown_audit
