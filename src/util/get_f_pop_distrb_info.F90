SUBROUTINE GET_F_HALO_MEMBER_VAR_ADDR

! use work_mod, only: haloMemberVarAddr
use domain,   only: POP_haloClinic

implicit none

  !  haloMemberVarAddr(1) = loc(POP_haloClinic%recvTask)
  !  haloMemberVarAddr(2) = loc(POP_haloClinic%sendTask)
  !  haloMemberVarAddr(3) = loc(POP_haloClinic%sizeSend)
  !  haloMemberVarAddr(4) = loc(POP_haloClinic%sizeRecv)
  !  haloMemberVarAddr(5) = loc(POP_haloClinic%srcLocalAddr)
  !  haloMemberVarAddr(6) = loc(POP_haloClinic%dstLocalAddr)
  !  haloMemberVarAddr(7) = loc(POP_haloClinic%sendAddr)
  !  haloMemberVarAddr(8) = loc(POP_haloClinic%recvAddr)

END SUBROUTINE

SUBROUTINE get_f_pop_distrb_info(numProcs, communicator, numLocalBlocks)

use LICOM_Error_mod
use precision_mod
use POP_DistributionMod
use POP_BlocksMod, only: POP_numBlocks
use domain, only: POP_distrbClinic, &
    POP_distrbClinic_blockLocation, &
    POP_distrbClinic_blockLocalID,  &
    POP_distrbClinic_blockGlobalID

implicit none
  integer (i4), intent(out) ::   &
      numProcs       ,&! number of processors in this dist
      communicator   ,&! communicator to use in this dist
      numLocalBlocks   ! number of blocks distributed to this
                       ! local processor

  integer (i4) :: errorCode ! returned error code
  call POP_DistributionGet(POP_distrbClinic, errorCode,      &
            numProcs       = numProcs,                       &
            communicator   = communicator,                   &
            numLocalBlocks = numLocalBlocks,                 &
            blockLocation  = POP_distrbClinic_blockLocation, &
            blockLocalID   = POP_distrbClinic_blockLocalID,  &
            blockGlobalID  = POP_distrbClinic_blockGlobalID)
END SUBROUTINE get_f_pop_distrb_info