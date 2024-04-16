SUBROUTINE fortran_mpi_barrier

use msg_mod, only:mpi_comm_ocn 
use param_mod, only: ierr

   IMPLICIT NONE

   call mpi_barrier(mpi_comm_ocn, ierr)

END SUBROUTINE fortran_mpi_barrier
