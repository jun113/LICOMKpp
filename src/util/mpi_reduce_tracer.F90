subroutine mpi_reduce_tracer(e1, e2)

use msg_mod, only: mpi_comm_ocn
use param_mod, only: ierr
use precision_mod, only: r8, mpi_pr, mpi_sum

implicit none
  real (r8) :: e1, e2

  call mpi_reduce(e1, e2, 1, mpi_pr, mpi_sum, 0, mpi_comm_ocn, ierr)

end subroutine mpi_reduce_tracer