SUBROUTINE get_block_info
use blocks
use domain

IMPLICIT NONE
   INTEGER      :: iblock
   type (block) :: this_block ! block information for current block
	 iblock = 1
   this_block = get_block(blocks_clinic(iblock),iblock)

	 ib = this_block%ib
	 ie = this_block%ie
	 jb = this_block%jb
	 je = this_block%je
END SUBROUTINE get_block_info

subroutine pop_halo_create_from_fortran_get(communicator, numMsgSend, numMsgRecv, numLocalCopies, &
                                             C_bufSendTask, C_bufRecvTask, &
                                             C_bufSizeSend, C_bufSizeRecv, &
                                             C_bufSrcLocalAddr, C_bufDstLocalAddr, &
                                             C_bufSendAddr, C_bufRecvAddr)
use pop_halomod
use precision_mod
use domain
implicit none
   integer (i4), intent(out) :: &
	 	 communicator, &
		 numMsgSend,   &
		 numMsgRecv,   &
		 numLocalCopies

   integer (i4), dimension(POP_haloClinic%numMsgSend), intent(out) :: &
      C_bufSendTask
   integer (i4), dimension(POP_haloClinic%numMsgRecv), intent(out) :: &
      C_bufRecvTask
   integer (i4), dimension(POP_haloClinic%numMsgSend), intent(out) :: &
      C_bufSizeSend
   integer (i4), dimension(POP_haloClinic%numMsgRecv), intent(out) :: &
      C_bufSizeRecv
   integer (i4), dimension(POP_haloClinic%numLocalCopies * 3), intent(out) :: &
      C_bufSrcLocalAddr
   integer (i4), dimension(POP_haloClinic%numLocalCopies * 3), intent(out) :: &
      C_bufDstLocalAddr
   integer (i4), dimension(POP_haloClinic%numMsgSend * 3 * bufSizeSend), intent(out) :: &
      C_bufSendAddr
   integer (i4), dimension(POP_haloClinic%numMsgRecv * 3 * bufSizeRecv), intent(out) :: &
      C_bufRecvAddr

   integer (i4) :: i, j, k ! dummy loop indices

	 communicator   = POP_haloClinic%communicator
	 numMsgSend     = POP_haloClinic%numMsgSend
	 numMsgRecv     = POP_haloClinic%numMsgRecv
	 numLocalCopies = POP_haloClinic%numLocalCopies

   ! if (.not. allocated(POPHalo_bufSendTask)) then
   !    allocate(POPHalo_bufSendTask(POP_haloClinic%numMsgSend))
	!  endif
   ! if (.not. allocated(POPHalo_bufRecvTask)) then
   !    allocate(POPHalo_bufRecvTask(POP_haloClinic%numMsgRecv))
	!  endif
   ! if (.not. allocated(POPHalo_bufSizeSend)) then
   !    allocate(POPHalo_bufSizeSend(POP_haloClinic%numMsgSend))
	!  endif
   ! if (.not. allocated(POPHalo_bufSizeRecv)) then
   !    allocate(POPHalo_bufSizeRecv(POP_haloClinic%numMsgRecv))
	!  endif
   ! if (.not. allocated(POPHalo_bufSrcLocalAddr)) then
   !    allocate(POPHalo_bufSrcLocalAddr(POP_haloClinic%numLocalCopies * 3))
	!  endif
   ! if (.not. allocated(POPHalo_bufDstLocalAddr)) then
   !    allocate(POPHalo_bufDstLocalAddr(POP_haloClinic%numLocalCopies * 3))
	!  endif
   ! if (.not. allocated(POPHalo_bufSendAddr)) then
   !    allocate(POPHalo_bufSendAddr(POP_haloClinic%numMsgSend * 3 * bufSizeSend))
	!  endif
   ! if (.not. allocated(POPHalo_bufRecvAddr)) then
   !    allocate(POPHalo_bufRecvAddr(POP_haloClinic%numMsgRecv * 3 * bufSizeRecv))
	!  endif

   do i=1, pop_haloclinic%numMsgSend
      C_bufSendTask(i) = POP_haloClinic%sendTask(i)
  	  C_bufSizeSend(i) = POP_haloClinic%sizeSend(i)
   enddo
   do i=1, pop_haloclinic%numMsgRecv
   	 C_bufRecvTask(i) = POP_haloClinic%recvTask(i)
   	 C_bufSizeRecv(i) = POP_haloClinic%sizeRecv(i)
   enddo
   do j=1, POP_haloClinic%numLocalCopies
   	 do i=1,3
   	 	 C_bufSrcLocalAddr((j-1)*3 + i) = POP_haloClinic%srcLocalAddr(i, j)
  	 	 C_bufDstLocalAddr((j-1)*3 + i) = POP_haloClinic%dstLocalAddr(i, j)
   	 enddo
   enddo
   do k=1, POP_haloClinic%numMsgSend
   	 do j=1, bufSizeSend
   	 	 do i=1, 3
   		 	 C_bufSendAddr((k-1) * (bufSizeSend * 3) + (j-1) * 3 + i) = POP_haloClinic%sendAddr(i, j, k)
   	 	 enddo
   	 enddo
   enddo
   do k=1, POP_haloClinic%numMsgRecv
   	 do j=1, bufSizeRecv
   	 	 do i=1, 3
   		 	 C_bufRecvAddr((k-1) * (bufSizeRecv * 3) + (j-1) * 3 + i) = POP_haloClinic%recvAddr(i, j, k)
   	 	 enddo
   	 enddo
   enddo

  ! write(*,*),"f pop", mytid, POP_haloClinic%communicator, POP_haloClinic%numMsgSend,pop_haloclinic%nummsgrecv,POP_haloClinic%numLocalCopies

end subroutine pop_halo_create_from_fortran_get

subroutine pop_halo_create_from_fortran_get_1(communicator, numMsgSend, numMsgRecv, numLocalCopies)
use pop_halomod
use precision_mod
use domain
implicit none
   integer (i4), intent(out) :: &
	 	 communicator, &
		 numMsgSend,   &
		 numMsgRecv,   &
		 numLocalCopies

   integer (i4) :: i, j, k ! dummy loop indices

	 communicator   = POP_haloClinic%communicator
	 numMsgSend     = POP_haloClinic%numMsgSend
	 numMsgRecv     = POP_haloClinic%numMsgRecv
	 numLocalCopies = POP_haloClinic%numLocalCopies

end subroutine pop_halo_create_from_fortran_get_1

subroutine pop_halo_create_from_fortran_get_2(C_bufSendTask, C_bufRecvTask, &
                                             C_bufSizeSend, C_bufSizeRecv, &
                                             C_bufSrcLocalAddr, C_bufDstLocalAddr, &
                                             C_bufSendAddr, C_bufRecvAddr)
use pop_halomod
use precision_mod
use domain
implicit none

   integer (i4), dimension(POP_haloClinic%numMsgSend), intent(out) :: &
      C_bufSendTask
   integer (i4), dimension(POP_haloClinic%numMsgRecv), intent(out) :: &
      C_bufRecvTask
   integer (i4), dimension(POP_haloClinic%numMsgSend), intent(out) :: &
      C_bufSizeSend
   integer (i4), dimension(POP_haloClinic%numMsgRecv), intent(out) :: &
      C_bufSizeRecv
   integer (i4), dimension(POP_haloClinic%numLocalCopies * 3), intent(out) :: &
      C_bufSrcLocalAddr
   integer (i4), dimension(POP_haloClinic%numLocalCopies * 3), intent(out) :: &
      C_bufDstLocalAddr
   integer (i4), dimension(POP_haloClinic%numMsgSend * 3 * bufSizeSend), intent(out) :: &
      C_bufSendAddr
   integer (i4), dimension(POP_haloClinic%numMsgRecv * 3 * bufSizeRecv), intent(out) :: &
      C_bufRecvAddr

   integer (i4) :: i, j, k ! dummy loop indices

   ! write(*,*), "F", loc(C_bufSendTask), loc(C_bufRecvTask), loc(C_bufSizeSend), loc(C_bufSizeRecv), loc(C_bufSrcLocalAddr), loc(C_bufDstLocalAddr), loc(C_bufSendAddr), loc(C_bufRecvAddr)

   do i=1, pop_haloclinic%numMsgSend
      C_bufSendTask(i) = POP_haloClinic%sendTask(i)
  	  C_bufSizeSend(i) = POP_haloClinic%sizeSend(i)
   enddo
   do i=1, pop_haloclinic%numMsgRecv
   	 C_bufRecvTask(i) = POP_haloClinic%recvTask(i)
   	 C_bufSizeRecv(i) = POP_haloClinic%sizeRecv(i)
   enddo
   do j=1, POP_haloClinic%numLocalCopies
   	 do i=1,3
   	 	 C_bufSrcLocalAddr((j-1)*3 + i) = POP_haloClinic%srcLocalAddr(i, j)
  	 	    C_bufDstLocalAddr((j-1)*3 + i)    = POP_haloClinic%dstLocalAddr(i, j)
   	 enddo
   enddo
   do k=1, POP_haloClinic%numMsgSend
   	 do j=1, bufSizeSend
   	 	 do i=1, 3
   		 	 C_bufSendAddr((k-1) * (bufSizeSend * 3) + (j-1) * 3 + i) = POP_haloClinic%sendAddr(i, j, k)
   	 	 enddo
   	 enddo
   enddo
   do k=1, POP_haloClinic%numMsgRecv
   	 do j=1, bufSizeRecv
   	 	 do i=1, 3
   		 	 C_bufRecvAddr((k-1) * (bufSizeRecv * 3) + (j-1) * 3 + i) = POP_haloClinic%recvAddr(i, j, k)
   	 	 enddo
   	 enddo
   enddo

  ! write(*,*),"f pop", mytid, POP_haloClinic%communicator, POP_haloClinic%numMsgSend,pop_haloclinic%nummsgrecv,POP_haloClinic%numLocalCopies

end subroutine pop_halo_create_from_fortran_get_2


! subroutine pop_halo_create_from_fortran_free_buf()
! use pop_halomod
! implicit none
!    if (allocated(POPHalo_bufSendTask)) then
!       deallocate(POPHalo_bufSendTask)
! 	 endif
!    if (allocated(POPHalo_bufRecvTask)) then
!       deallocate(POPHalo_bufRecvTask)
! 	 endif
!    if (allocated(POPHalo_bufSizeSend)) then
!       deallocate(POPHalo_bufSizeSend)
! 	 endif
!    if (allocated(POPHalo_bufSizeRecv)) then
!       deallocate(POPHalo_bufSizeRecv)
! 	 endif
!    if (allocated(POPHalo_bufSrcLocalAddr)) then
!       deallocate(POPHalo_bufSrcLocalAddr)
! 	 endif
!    if (allocated(POPHalo_bufDstLocalAddr)) then
!       deallocate(POPHalo_bufDstLocalAddr)
! 	 endif
!    if (allocated(POPHalo_bufSendAddr)) then
!       deallocate(POPHalo_bufSendAddr)
! 	 endif
!    if (allocated(POPHalo_bufRecvAddr)) then
!       deallocate(POPHalo_bufRecvAddr)
! 	 endif
! end subroutine pop_halo_create_from_fortran_free_buf