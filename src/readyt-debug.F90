!  CVS: $Id: readyt.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE READYT_debug(wkk1,wkk2,wkk3,wkk4,wkk5,wkk6,wkk7,wkk8,wkk9,wkk10,wkk11,wkk12,wkk13,wkk14)
!     =================
!     PREPARATION OF BAROTROPIC AND BAROCLINIC INTEGRATION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use dyn_mod
use work_mod
use pmix_mod
use msg_mod
use forc_mod, only: psa,USTAR,BUOYTUR, BUOYSOL,NSWV,SWV,su,sv
use domain
use grid
use blocks
use constant_mod
use operators

!
      IMPLICIT NONE
      real(r8), dimension(imt,jmt,km)::wkk1,wkk2,wkk3,wkk4,wkk7,wkk8,wkk9,wkk10,wkk11,wkk13,wkk14
      real(r8), dimension(imt,jmt,km,ntra)::wkk5
      real(r8), dimension(imt,jmt,kmm1)::wkk6,wkk12
      REAL(r8)   :: ABCD,TUP,SUP,TLO,SLO,RHOUP,RHOLO,ek0
      REAL(r8)   :: DENS, zzz1,zzz2,epsln
      real(r8),dimension(:,:,:,:),allocatable:: alpha, beta, pp, ppa, ppb, ppc
      real(r8),dimension(:,:,:),allocatable:: work1, work2, adv_tt
      integer :: iblock
      EXTERNAL DENS
      type (block) :: this_block          ! block information for current block
 
      epsln = 1.0D-25 
      allocate ( alpha(imt,jmt,km,max_blocks_clinic), beta(imt,jmt,km,max_blocks_clinic))
      allocate ( pp(imt,jmt,km,max_blocks_clinic), ppa(imt,jmt,km,max_blocks_clinic))
      allocate ( ppb(imt,jmt,km,max_blocks_clinic), ppc(imt,jmt,km,max_blocks_clinic))
      allocate ( work1(imt,jmt,max_blocks_clinic), work2(imt,jmt,max_blocks_clinic),adv_tt(imt,jmt,km))

      allocate(dlu(imt,jmt,km,max_blocks_clinic),dlv(imt,jmt,km,max_blocks_clinic),gg(imt,jmt,km,max_blocks_clinic))
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
    do iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UTL (I,J,K,IBLOCK)= wkk1 (I,J,K)
               VTL (I,J,K,IBLOCK)= wkk2 (I,J,K)
               UTF (I,J,K,IBLOCK)= wkk3 (I,J,K)
               VTF (I,J,K,IBLOCK)= wkk4 (I,J,K)
               pdensity (I,J,K,IBLOCK)= wkk7 (I,J,K)
               gg (I,J,K,IBLOCK)= wkk8 (I,J,K)
               pp (I,J,K,IBLOCK)= wkk9 (I,J,K)
               alpha (I,J,K,IBLOCK)= wkk10 (I,J,K)
               beta (I,J,K,IBLOCK)= wkk11 (I,J,K)
               dlu (I,J,K,IBLOCK)= wkk13 (I,J,K)
               dlv (I,J,K,IBLOCK)= wkk14 (I,J,K)
            END DO
         END DO
      END DO
    end do
    do iblock = 1, ntra
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               akt (I,J,K,IBLOCK,1)= wkk5 (I,J,K,iblock)
            END DO
         END DO
      END DO
    end do
      allocate(rict(imt,jmt,kmm1,max_blocks_clinic))
    do iblock = 1, nblocks_clinic
      DO K = 1,KMm1
         DO J = JST,JET
            DO I = 1,IMT
               rict (I,J,K,IBLOCK)= wkk6 (I,J,K)
               ricdt (I,J,K,IBLOCK)= wkk12 (I,J,K)
            END DO
         END DO
      END DO
    end do
    do iblock = 1, nblocks_clinic
         DO J = JST,JET
            DO I = 1,IMT
               wka (I,J,1,IBLOCK)= su (I,J,iblock)
               wka (I,J,2,IBLOCK)= sv (I,J,iblock)
            END DO
         END DO
    end do
 
      allocate(rit(imt,jmt,kmm1,max_blocks_clinic),ric(imt,jmt,kmm1,max_blocks_clinic))
      allocate(rict_replace(imt,jmt,kmm1,max_blocks_clinic))

    rit   = 0.0_r8
    ric   = 0.0_r8
   !  rict  = 0.0_r8
    ricdtTmS  = 0.0_r8
   !  ricdt = 0.0_r8
   rict_replace = rict

      deallocate ( alpha, beta)
      deallocate ( pp, ppa, ppb,ppc)
      deallocate ( work1, work2, adv_tt)
  call mpi_barrier(mpi_comm_ocn,ierr)
!
      if (ist == 0 ) then
         ax = c0
         ay = c0
         az = c0
      end if
!
      RETURN
      END SUBROUTINE READYT_debug
 
 
