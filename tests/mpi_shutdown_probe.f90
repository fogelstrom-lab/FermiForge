program mpi_shutdown_probe
  use mpi_f08
  use mpi_shutdown_audit, only: finalize_with_audit
  use iso_c_binding, only: c_int
  implicit none
  interface
    function c_sleep(seconds) bind(C,name="sleep") result(left)
      import c_int
      integer(c_int), value :: seconds
      integer(c_int) :: left
    end function
  end interface
  integer :: ierr, rank, unit, left
  character(1024) :: prefix
  character(32) :: suffix
  call get_command_argument(1,prefix)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  write(suffix,'(".ready.",i0)') rank
  open(newunit=unit,file=trim(prefix)//trim(suffix),status='replace')
  write(unit,*) 'ready'
  close(unit)
  left=c_sleep(3_c_int)
  call finalize_with_audit(trim(prefix),ierr)
  if (ierr /= MPI_SUCCESS) error stop 'finalize failed'
end program
