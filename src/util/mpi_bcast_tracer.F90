subroutine mpi_bcast_tracer(err_norm)

use msg_mod, only: mpi_comm_ocn
use param_mod, only: ierr
use precision_mod, only: r8, mpi_pr

implicit none

   real (r8) :: err_norm

   call mpi_bcast(err_norm, 1, mpi_pr, 0, mpi_comm_ocn, ierr)

end subroutine