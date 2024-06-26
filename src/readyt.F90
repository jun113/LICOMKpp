!  CVS: $Id: readyt.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE READYT
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
use forc_mod, only: psa,USTAR,BUOYTUR, BUOYSOL,NSWV,SWV
use domain
use grid
use blocks
use constant_mod
use operators

!
      IMPLICIT NONE
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

!     if (ist == 0 .and. iday ==2) then
!       do iblock = 1, nblocks_clinic
!       do k=1,km
!       do j=1 , jmt
!       do i = 1,imt
!          at(i,j,k,1,iblock) = to(k)
!          at(i,j,k,2,iblock) = so(k)
!          atb(i,j,k,1,iblock) = to(k)
!          atb(i,j,k,2,iblock) = so(k)
!       end do
!       end do
!       end do
!       end do
!     end if
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    do iblock = 1, nblocks_clinic
      DO J = JST,JET
         DO I = 1,IMT
            H0L (I,J,IBLOCK)= H0F (I,J,IBLOCK)
            H0F (I,J,IBLOCK)= H0 (I,J,IBLOCK)
         END DO
      END DO
    end do
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
    do iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UTL (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK)
               VTL (I,J,K,IBLOCK)= VTF (I,J,K,IBLOCK)
               UTF (I,J,K,IBLOCK)= U (I,J,K,IBLOCK)
               VTF (I,J,K,IBLOCK)= V (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
    end do
 
!     --------------------------------------------------------------
!     PREPARING FOR CALCULATING RICHARDSON NUMBER
!     --------------------------------------------------------------
      allocate(dlu(imt,jmt,km,max_blocks_clinic),dlv(imt,jmt,km,max_blocks_clinic),gg(imt,jmt,km,max_blocks_clinic))

!      callocate in pmix_mod ric rict rict_replace wjl 20210826
!      allocate(rit(imt,jmt,kmm1,max_blocks_clinic))
      allocate(rit(imt,jmt,kmm1,max_blocks_clinic),ric(imt,jmt,kmm1,max_blocks_clinic))
      allocate(rict(imt,jmt,kmm1,max_blocks_clinic))
      allocate(rict_replace(imt,jmt,kmm1,max_blocks_clinic))

!      allocate(rict_ref(imt,jmt,max_blocks_clinic))
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,K)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
!lhl20160425          call tgrid_to_ugrid(dlu(:,:,k,iblock), at(:,:,k,1,iblock),iblock)
!lhl20160425          call tgrid_to_ugrid(dlv(:,:,k,iblock), at(:,:,k,2,iblock),iblock)
          call tgrid_to_ugrid(dlu(:,:,k,iblock), atb(:,:,k,1,iblock),iblock)
          call tgrid_to_ugrid(dlv(:,:,k,iblock), atb(:,:,k,2,iblock),iblock)
      END DO
   END DO

    rit   = 0.0_r8
    ric   = 0.0_r8
    rict  = 0.0_r8
    ricdtTmS  = 0.0_r8
    ricdt = 0.0_r8
    akt   = 0.0_r8
!
! ric locate at integer level
!

!$OMP PARALLEL DO PRIVATE(IBLOCK,TUP,SUP,TLO,SLO,RHOUP,RHOLO)
   DO IBLOCK = 1,NBLOCKS_CLINIC
      DO K = 1,KMM1
         DO J = 1, JMT-1
            DO I = 2,IMT
               TUP = DLU (I,J,K,IBLOCK) - TO (K +1)
               SUP = DLV (I,J,K,IBLOCK) - SO (K +1)
               TLO = DLU (I,J,K +1,IBLOCK) - TO (k +1)
               SLO = DLV (I,J,K +1,IBLOCK) - SO (K +1)
               RHOUP = DENS (TUP, SUP, K +1)
               RHOLO = DENS (TLO, SLO, K +1)
               ric (I,J,K,IBLOCK) = VIV (I,J,K +1,IBLOCK)* OD0* G * (RHOLO - RHOUP)*ODZT(K+1)
            END DO
         END DO
      END DO
    END DO
!
! calculate the potential density using BC equation on U-grid
! note the referrence of density was added!
! if the salinity should multiple 1000??????
      CALL DENSITY
!
! ric locate at integer level from 2
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,TUP,SUP,TLO,SLO,RHOUP,RHOLO)
   DO IBLOCK = 1,NBLOCKS_CLINIC
      DO K = 1,KMM1
         DO J = 1,JMT
            DO I = 1, IMT
!lhl20160425               TUP = AT (I,J,K,1,IBLOCK) - TO (K +1)
!lhl20160425               SUP = AT (I,J,K,2,IBLOCK) - SO (K +1)
!lhl20160425               TLO = AT (I,J,K +1,1,IBLOCK) - TO (k +1)
!lhl20160425               SLO = AT (I,J,K +1,2,IBLOCK) - SO (K +1)
               TUP = ATB (I,J,K,1,IBLOCK) - TO (K +1)
               SUP = ATB (I,J,K,2,IBLOCK) - SO (K +1)
               TLO = ATB (I,J,K +1,1,IBLOCK) - TO (k +1)
               SLO = ATB (I,J,K +1,2,IBLOCK) - SO (K +1)
               RHOUP = DENS (TUP, SUP, K +1)
               RHOLO = DENS (TLO, SLO, K +1)
               rict (I,J,K,IBLOCK) = VIT (I,J,K +1,IBLOCK)* OD0* G * (RHOLO - RHOUP)*ODZT(K+1) !LPF20160730
              ! rict (I,J,K,IBLOCK) = VIT (I,J,K +1,IBLOCK)* OD0* G * (PDENSITY (I,J,K+1,IBLOCK) - PDENSITY (I,J,K,IBLOCK))*ODZT(K+1)
            END DO
         END DO
      END DO
   END DO

!LPF20160516
!Move to readyc due to caculation of MLD 
!LPF20160516
   rict_replace = rict

!
!     --------------------------------------------------------------
!     COMPUTING DENSITY AND BAROCLINIC PRESSURE
!     --------------------------------------------------------------
 
! calculate the potential density using BC equation on U-grid
! note the referrence of density was added!
! if the salinity should multiple 1000??????
!      CALL DENSITY !move to the upper due to using the density

!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               GG (I,J,K,IBLOCK)=-OD0*G*PDENSITY (I,J,K,IBLOCK)*VIT(I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO

!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1, JMT
         DO I = 1,IMT
! calculate pressure at level 1
            PP (I,J,1,IBLOCK)= GG (I,J,1,IBLOCK)*0.5D0* DZP (1)* VIT (I,J,1,IBLOCK)
!
            PPA(I,J,1,IBLOCK)= PSA (I,J,IBLOCK)* VIT (I,J,1,IBLOCK)
            PPB(I,J,1,IBLOCK)= AT(I,J,1,1,IBLOCK)*VIT(I,J,1,IBLOCK)
            PPC(I,J,1,IBLOCK)= AT(I,J,1,2,IBLOCK)*VIT(I,J,1,IBLOCK)
         END DO
      END DO
   END DO


!$OMP PARALLEL DO PRIVATE (IBLOCK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 2,KM
         DO J = 1, JMT
            DO I = 1,IMT
!lhl1204
! calculate pressure at T-grid
               PP (I,J,K,IBLOCK)= VIT (I,J,K,IBLOCK)* (PP(I,J,K -1,IBLOCK) +0.5D0* &
               (GG (I,J,K,IBLOCK)* DZP (K) + GG (I,J,K -1,IBLOCK)* DZP (K -1)))
               PPA(I,J,K,IBLOCK)= VIT (I,J,K,IBLOCK)* (PPA(I,J,K -1,IBLOCK)+GG (I,J,K -1,IBLOCK)*DZP (K -1))
               PPB(I,J,K,IBLOCK)= VIT (I,J,K,IBLOCK)* (AT (I,J,K-1,1,IBLOCK)-  &
                                  (AT (I,J,K-1,1,IBLOCK)-AT (I,J,K,1,IBLOCK))*DZP(K -1)/(DZP(K-1)+DZP(K)))
               PPC(I,J,K,IBLOCK)= VIT (I,J,K,IBLOCK)* (AT (I,J,K-1,2,IBLOCK)-  &
                                  (AT (I,J,K-1,2,IBLOCK)-AT (I,J,K,2,IBLOCK))*DZP(K -1)/(DZP(K-1)+DZP(K)))
            END DO
         END DO
      END DO
   END DO
!
!lhl0711#endif
!lhl1204
! calculate the thermal expansion and the salinity contraction at T-grid
!  PP in negative in model
!  in PP a OD0 is mutilped and should be divided
!  PP is in Pa to tranform to db by divided 10000.
!
!      CALL THERMAL(AT(1,1,1,1),AT(1,1,1,2),PP,ALPHA,BETA,VIT)
      CALL THERMAL(PPB,PPC,PPA,ALPHA,BETA,VIT)
!
! calculate the surface buoyancy fluxes at T-grid

!M
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1, JMT
         DO I = 1,IMT
         BUOYTUR(I,J,IBLOCK)= VIT(I,J,1,IBLOCK)*NSWV(I,J,IBLOCK)*G*ALPHA(I,J,1,IBLOCK)*OD0CP
         BUOYSOL(I,J,IBLOCK)= VIT(I,J,1,IBLOCK)* SWV(I,J,IBLOCK)*G*ALPHA(I,J,1,IBLOCK)*OD0CP
         END DO
      END DO
   END DO
!
!
! calculate the T component minus S component of B-V frequency at T-grid
!

!M
!$OMP PARALLEL DO PRIVATE (iblock)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KMM1
         DO J = 1, JMT
            DO I = 2, IMT-1
               ricdtTmS(I,J,K,iblock)=VIT(I,J,K+1,iblock)*G*((AT(I,J,K,1,iblock)-AT(I,J,K+1,1,iblock)) &
                               *ALPHA(I,J,K+1,iblock)+1000.D0*(AT(I,J,K,2,iblock)-AT(I,J,K+1,2,iblock)) &
                               *BETA(I,J,K+1,iblock))*ODZT(K+1)
!changed the form by LPF20160728
               ricdt(I,J,K,iblock)=VIT(I,J,K+1,iblock)/((AT(I,J,K,1,iblock)-AT(I,J,K+1,1,iblock)+epsln)* & 
                                   ALPHA(I,J,K+1,iblock))*1000.D0*((AT(I,J,K,2,iblock)-AT(I,J,K+1,2,iblock))*  &
                                   BETA(I,J,K+1,iblock))
               !ricdt(I,J,K,iblock)=(AT(I,J,K,2,iblock)-AT(I,J,K+1,2,iblock))*BETA(I,J,K+1,iblock)*1000.D0 &
               !                    /((AT(I,J,K,1,iblock)-AT(I,J,K+1,1,iblock)+epsln)* & 
               !                    ALPHA(I,J,K+1,iblock))*VIT(I,J,K+1,iblock)
!changed the form by LPF20160728
#if ( defined CANUTOMIXOUT )
               alpha_canuto(I,J,K,iblock)=ALPHA(I,J,K,iblock)     !yuzp-2016/12/4
               beta_canuto(I,J,K,iblock)=BETA(I,J,K,iblock)     !yuzp-2016/12/4
#endif
            END DO
         END DO
      END DO
   END DO

!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               GG (I,J,K,IBLOCK)=-OD0*G*(PDENSITY(I,J,K,IBLOCK)-PO(K)-1000.0D0)*VIT(I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
!lhl1204
 
!     --------------------------------------------------------------
!     COMPUTING HORIZONTAL GRADIENT OF BAROCLINIC PRESSURE
!     --------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (IBLOCK, this_block)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call grad(k, dlu(:,:,k,iblock), dlv(:,:,k,iblock), PP(:,:,k,iblock), this_block)
      END DO
   END DO
 
!     --------------------------------------------------------------
!     COMPUTING VERTICALLY INTEGRATED GRADIENT OF BAROCLINIC PRESSURE
!     --------------------------------------------------------------
 
      CALL VINTEG (DLU,PXB)
      CALL VINTEG (DLV,PYB)

 
!     --------------------------------------------------------------
!     COMPUTING WGP, WHX, WHY, PAX & PAY
!     --------------------------------------------------------------
  dlu(:,:,1,:) = 0.0D0
  dlu(:,:,2,:) = 0.0D0
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,ABCD)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
        DO J = 1, JMT
         DO I = 1,IMT
               ABCD = GG (I,J,K,IBLOCK)* OHBT (I,J,IBLOCK)* DZP (K)
               DLU (I,J,1,IBLOCK)= DLU (I,J,1,IBLOCK) + ABCD
               DLU (I,J,2,IBLOCK)= DLU (I,J,2,IBLOCK) + ABCD * ZKT (K)
            END DO
         END DO
      END DO
  END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1, JMT
         DO I = 1,IMT
            DLV (I,J,1,IBLOCK)= (DLU (I,J,1,IBLOCK) + DLU (I,J,2,IBLOCK)* OHBT (I,J,IBLOCK))/ G
            DLV (I,J,2,IBLOCK)= DLU (I,J,2,IBLOCK)* OHBT (I,J,IBLOCK)* OHBT (I,J,IBLOCK)
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      call tgrid_to_ugrid(wgp(:,:,iblock), dlv(:,:,1,iblock),iblock)
      call tgrid_to_ugrid(work(:,:,iblock), dlv(:,:,2,iblock),iblock)
      whx(:,:,iblock) = hbx(:,:,iblock)*work(:,:,iblock)*viv(:,:,1,iblock)
      why(:,:,iblock) = hby(:,:,iblock)*work(:,:,iblock)*viv(:,:,1,iblock)
  END DO

!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      call grad(1, work1(:,:,iblock), work2(:,:,iblock), psa , this_block)
   END DO
   pay = -OD0*work2
   pax = -OD0*work1
!!
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
      END SUBROUTINE READYT
 
 
