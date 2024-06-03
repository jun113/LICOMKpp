!     =================
      SUBROUTINE READYC_debug(wkk1,wkk2,wkk3,wkk4,wkk5,wkk6,wkk7,wkk8,wkk9,wkk10,wkk11,wkk12,wkk13)
!     =================
!     PREPARATION OF BAROTROPIC AND BAROCLINIC INTEGRATION
 
!     ADVECTION + DIFFUSION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use tracer_mod
use pmix_mod
#if ( defined TIDEMIX )
use forc_mod, only: su,sv,USTAR,BUOYTUR, BUOYSOL,wave_dis
#else
use forc_mod, only: su,sv,USTAR,BUOYTUR, BUOYSOL
#endif
use domain
use grid
use blocks
use advection
use operators
use LICOM_Error_mod
#ifdef BIHAR
use hmix_del4
#else
use hmix_del2
#endif
use msg_mod 
use gather_scatter
use distribution
use constant_mod
use canuto_2010_mod
      IMPLICIT NONE
!      REAL(r8)  :: WKP (KMP1)
      REAL(r8),dimension(IMT,JMT,KM) :: wkk1,wkk2,wkk5,wkk6,wkk8,wkk9,wkk10,wkk11
      REAL(r8),dimension(IMT,JMT) :: wkk7,wkk12,wkk13
      REAL(r8),dimension(IMT,JMT,KMm1) :: wkk3,wkk4
      INTEGER   :: IWK,n2, iblock,kmb
      REAL(r8)  :: WK1 (KM) ,WK2 (KM), WK3 (KM),WK4(KM)
      REAL(r8)  :: WP1 (KM) ,WP2 (KM), WP3 (KM)
      REAL(r8)  :: WP4 (KMM1) ,WP5 (KMM1), WP6 (KMM1)
      REAL(r8)  :: WP7 (KMM1) ,WP8(KM)
      REAL(r8)  :: WP9 ,WP10, WP11
      REAL(r8),dimension(IMT,JMT,KM,MAX_BLOCKS_CLINIC) :: WP12,WP13,riv1,riv2
      REAL(r8)  :: epsln,RKV,RKV1,ek0
      REAL(r8)  :: adv_x1,adv_x2,adv_x,adv_y1,adv_y2,adv_z,diff_u1,diff_u2,diff_v1,diff_v2
      REAL(r8)  :: dlux,dlvx,dluy,dlvy,dluz,dlvz,adv_z1,adv_z2,adv_z3,adv_z4
      real(r8)  :: hdvk(imt,jmt), hduk(imt,jmt), adv_uu(imt,jmt,km), adv_vv(imt,jmt,km)
      real(r8)  :: tq,sq,ak_tide_mixing(km)
#ifdef BCKMEX
      real(r8) :: diff_back(imt,jmt,max_blocks_clinic),&
                  diff_back_sh(imt,jmt,max_blocks_clinic),&
                  diff_back_nh(imt,jmt,max_blocks_clinic)
#endif
      REAL(r8)    :: DENS
      EXTERNAL DENS

!
!     real (r8) :: ttt(imt_global, jmt_global)
      type (block) :: this_block          ! block information for current block
       
#ifdef CANUTO2010
!#if (defined CANUTO2010)
      REAL(r8)  :: AIDIF
#endif

!YU 
      real (r8):: akt_back(2*km),aks_back(2*km),akm_back(2*km)
      allocate(riu(imt,jmt,0:km,max_blocks_clinic),stat=ierr)
         if (.not. allocated(rict_ref)) allocate(rict_ref(imt,jmt,max_blocks_clinic))

      allocate(dlub(imt,jmt,max_blocks_clinic),dlvb(imt,jmt,max_blocks_clinic),stat=ierr)
!Yu   if (ierr /= 0) then
!Yu      write(6,*)'allocation error---riu'
!Yu      stop
!Yu   end if
 
      epsln = 1.0D-25
! !$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
!    DO IBLOCK = 1, NBLOCKS_CLINIC
!       DO J = JST,JET
!          DO I = 1,IMT
!             H0BL (I,J,IBLOCK)= H0BF (I,J,IBLOCK)
!             H0BF (I,J,IBLOCK)= H0 (I,J,IBLOCK)
!          END DO
!       END DO
!    END DO
 
!lhl0711
#if (defined CANUTO2010)
       AIDIF=0.5  
!      if (ISC/=0)  AIDIF=0.5  
#endif
!lhl0711
 
!---------------------------------------------------------------------
!     Calculating Richardson number riu (at U/V-point);
!---------------------------------------------------------------------
      s2t  = c0
      ridt = c0
      riu  = c0
      wp12 = c0
      wp13 = c0

    do iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               wp12 (I,J,K,IBLOCK)= wkk1 (I,J,K)
               wp13 (I,J,K,IBLOCK)= wkk2 (I,J,K)
               ws (I,J,K,IBLOCK)= wkk6 (I,J,K)
               wka (I,J,K,IBLOCK)= wkk5 (I,J,K)
               dlu (I,J,K,IBLOCK)= wkk10 (I,J,K)
               dlv (I,J,K,IBLOCK)= wkk11 (I,J,K)
            END DO
         END DO
      END DO
    end do
    do iblock = 1, nblocks_clinic
      DO K = 1,KMm1
         DO J = JST,JET
            DO I = 1,IMT
               rit (I,J,K,IBLOCK)= wkk3 (I,J,K)
               s2t (I,J,K,IBLOCK)= wkk4 (I,J,K)
            END DO
         END DO
      END DO
    end do
    do iblock = 1, nblocks_clinic
         DO J = JST,JET
            DO I = 1,IMT
               dlub (I,J,IBLOCK)= wkk12 (I,J)
               dlvb (I,J,IBLOCK)= wkk13 (I,J)
            END DO
         END DO
    end do

! #ifdef CANUTO2010
!
   !   amld = c0
   !   akmt = c0
   !   akmu = c0
!

! !$OMP PARALLEL DO PRIVATE (iblock,tq,sq,wp1,wp2,wp3,wp4,wp5,wp6,wp7,wp8,wk1,wk2,wk3,wk4)
!    do iblock = 1, nblocks_clinic
!       DO J = 2, JMT-1
!          DO I = 2, IMT-1

!         if (VIT(I,J,1,iblock).gt.0.5) then
! !
!          wp1=0.D0
!          wp2=0.D0
!          wp3=0.D0
!          wp4=0.D0
!          wp5=0.D0
!          wp6=0.D0
!          wp7=0.D0
!          wp8=0.D0
! !
!          do k=1,km - 1
!             wp8(k)=-vit(i,j,k+1,iblock)*ZKP(k+1)
!          end do
! !
!          DO K = 1,KMT(i,j,iblock)-1
!                wp1(K)= VIT (I,J,K+1,iblock)* (AT (I,J,K,1,iblock)-(AT (I,J,K,1,iblock)-AT (I,J,K+1,1,iblock))* &
!                                DZP (K)/(DZP(K)+DZP(K+1)))
!                wp2(K)= VIT (I,J,K+1,iblock)* (AT (I,J,K,2,iblock)-(AT (I,J,K,2,iblock)-AT (I,J,K+1,2,iblock))* &
!                                DZP (K)/(DZP(K)+DZP(K+1)))*1.0D3+35.
!                wp4(k)=vit(i,j,k+1,iblock)*RIT(i,j,k,iblock)
!                wp5(k)=vit(i,j,k+1,iblock)*RICDT(i,j,k,iblock)
!                wp6(k)=vit(i,j,k+1,iblock)*S2T(i,j,k,iblock)
!                wp7(k)=vit(i,j,k+1,iblock)*RICT(i,j,k,iblock)

!                if (mytid == 3) then
!                   ! write(*,*)k,j,i,wp1(k) - wkk5(i,j,k)
!                   ! if (at(i,j,k,1,iblock) - wkk5(i,j,k) /= 0.0) then
!                   ! if (wp1(k) - wkk5(i,j,k) /= 0.0) then
!                   ! write(*,*)k,j,i,at(i,j,k,1,iblock) - wkk5(i,j,k)
!                   ! write(*,*)k,j,i,wp1(k) - wkk5(i,j,k)
!                   ! endif
!                end if
!          END DO
! !
!          DO K = 1,KMT(i,j,iblock)
!                tq = at(i,j,k,1,iblock) - to(1)
!                sq = at(i,j,k,2,iblock) - so(1)
!                wp3(k) = dens(tq,sq,1)+1.0D3
!          end do
! !
! !YU  May 10th, 2017
!          where ( wp7 < 0.0 )  wp7 = 0.0_r8
! !YU  May 10th, 2017
!          call  canuto_2010_interface_debug(wk1,    wk2,     wk3,     wk4,     amld(i,j,iblock) ,&
!                                      wp1,    wp2,     wp3,     wp4,     wp5,              &
!                                      wp7,    wp6,     ulat(i,j,iblock)/degtorad ,     wp8,& 
!                                      kmt(i,j,iblock) ,i,j, iblock, isc,ustar(i,j,iblock),wkk5)

! #if ( defined TIDEMIX )
          
!           ak_tide_mixing(:)=0.0
!           !LPF2016Nov12 IF (1./OHBT(I,J,IBLOCK).GT.SHELF_CUTOFF) THEN
!           !IF (1./OHBT(I,J,IBLOCK).GT.SHELF_CUTOFF) THEN
!           DO K = 1,INT(KMT(I,J,IBLOCK))-1
!               ak_tide_mixing(K)=back_tidalmixing+mixing_ef*local_mixing_fraction* & 
!                                     WAVE_DIS(I,J,IBLOCK)*FZ_TIDE(I,J,K,IBLOCK) &
!                                    /(dmax1(rict(I,J,K,iblock),1.0d-8)*WP3(K))     !yuzp-2016/11/23--12/2
!               ak_tide_mixing(K)=DMIN1(ak_tide_mixing(K),max_tidalmixing)   !LPF2016Nov12
!               richardson(I,J,K,iblock)=rict(I,J,K,IBLOCK)     !yuzp-2016/11/13--12/2
!               fztidal(I,J,K,iblock)=FZ_TIDE(I,J,K,IBLOCK)     !yuzp-2016/11/13
!               wp3_tidal(I,J,K,iblock)=WP3(K)     !yuzp-2016/11/13
!           END DO
          
!           DO K = INT(KMT(I,J,IBLOCK))-2,1,-1     !yuzp-2016/11/23
!            ak_tide_mixing(K)=DMIN1(ak_tide_mixing(K),ak_tide_mixing(K+1))     !yuzp-2016/11/23
!           END DO     !yuzp-2016/11/23
! #endif
!          DO K = 1,KM
!                ! if((mytid == 3) .and. (wk1(k)-wkk5(i,j,k) /= 0.0)) then
!               !  if(mytid == 3) then
!                ! write(*,*) k,j,i,wk1(k)-wkk5(i,j,k)
!                ! endif
!          !   AKMT(I,J,K,iblock)=WK1(K)
!          !   AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+ WK2(K)
!          !   AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+ WK3(K)
! #if ( defined TIDEMIX )
!          !   AKMT(I,J,K,iblock)=AKMT(I,J,K,iblock)+  ak_tide_mixing(K)*5.0
!          !   AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+ak_tide_mixing(K)
!          !   AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+ak_tide_mixing(K)
! #endif 
!          !   akmt(i,j,k,iblock) = min(akmt(i,j,k,iblock), 8.0D-3*DZP(K)*DZP(K))
!          !   akt(i,j,k,1,iblock)= min(akt(i,j,k,1,iblock),8.0D-3*DZP(K)*DZP(K))
!          !   akt(i,j,k,2,iblock)= min(akt(i,j,k,2,iblock),8.0D-3*DZP(K)*DZP(K))
!          END DO
!         endif
!          END DO
!       END DO
!    END DO
! !  stop
! !!
! ! calculate the vertical mixing on U-grid
!    ! call fortran_mpi_barrier
!    ! stop

! !$OMP PARALLEL DO PRIVATE (IBLOCK)
!    ! do iblock = 1, nblocks_clinic
!    !    DO K = 1,KMM1
!    !       call tgrid_to_ugrid(akmu(:,:,k,iblock), akmt(:,:,k,iblock),iblock)
!    !       do j= 1,jmt
!    !       do i= 1,imt
!    !          akmu(i,j,k,iblock) = akmu(i,j,k,iblock)* viv(i,j,k+1,iblock)
!    !       end do
!    !       end do
!    !    END DO
!    ! end do
! !lhl241204
! #endif
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM: ZONAL COMPONENT
!---------------------------------------------------------------------
      ! CALL UPWELL_debug (U,V,H0,wkk7,wkk8,wkk9,wkk5)
      
      ! call fortran_mpi_barrier
      ! stop
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
   !   dlu = c0
   !   dlv = c0
   !   wka = c0

! !$OMP PARALLEL DO PRIVATE (IBLOCK,K)
!    do iblock = 1, nblocks_clinic
!        DO K = 1,KM
!           call tgrid_to_ugrid(wka(:,:,k,iblock),ws(:,:,k,iblock), iblock)
!        END DO
!    end do

!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERMS
!---------------------------------------------------------------------
 
! !$OMP PARALLEL DO PRIVATE (IBLOCK)
!     DO IBLOCK = 1, NBLOCKS_CLINIC
!         call advection_momentum(u(:,:,:,iblock),v(:,:,:,iblock),wka(:,:,:,iblock), &
!                                 dlu(:,:,:,iblock),dlv(:,:,:,iblock),iblock)
! !       dlu(:,:,:,IBLOCK)=adv_uu
! !       dlv(:,:,:,IBLOCK)=adv_vv
!    END DO
!
!---------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL VISCOSITY
!---------------------------------------------------------------------
 
! #if ( defined SMAG)
! !
!          CALL SMAG2 (K)
! !
! #if (defined SMAG_FZ )
! #else
! #endif
! #else
! #if (defined BIHAR)

! !$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,HDUK,HDVK)
!    DO IBLOCK = 1, NBLOCKS_CLINIC
!       this_block = get_block(blocks_clinic(iblock),iblock)
!       DO K = 1,KM
!           call hdiffu_del4(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
!           do j = 3, jmt-2
!           do i = 3, imt-2
!             dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
!             dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
!           end do
!           end do
!       END DO
!    END DO

! #else
! !$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,hduk,hdvk,K)
!    DO IBLOCK = 1, NBLOCKS_CLINIC
!       this_block = get_block(blocks_clinic(iblock),iblock)
!       DO K = 1,KM
!          call hdiffu_del2(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
!          DO J = 3, JMT-2
!             DO I = 3,IMT-2
!                dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
!                dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
!             END DO
!          END DO
!       END DO
!    END DO
 
! #endif
! #endif
 
 
!---------------------------------------------------------------------
!     VERTICAL INTEGRATION
!---------------------------------------------------------------------
      ! CALL VINTEG (DLU,DLUB)
      ! CALL VINTEG (DLV,DLVB)
!
! !$OMP PARALLEL DO PRIVATE (iblock,J,I,kmb) 
!       do iblock = 1, nblocks_clinic
!          DO J = 2, jmt-1
!             DO I = 2,imt-1
!                kmb = kmu(i,j,iblock)
!                SBCX(I,J,IBLOCK) = SU (I,J,IBLOCK)* OD0
!                SBCY(I,J,IBLOCK) = SV (I,J,IBLOCK)* OD0
!                BBCX(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+  &
!                                           VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))   &
!                           *(UP(I,J,kmb,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP (I,J,kmb,IBLOCK)*SAG)
!                BBCY(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+  &  
!                                           VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))&
!                           *(-SNLAT(I,J,IBLOCK)*UP(I,J,kmb,IBLOCK)*SAG+VP(I,J,kmb,IBLOCK)*CAG)
!               dlub(i,j,iblock) = dlub(i,j,iblock) +(SBCX(i,j,iblock)-BBCX(i,j,iblock))*ohbu(i,j,iblock)
!               dlvb(i,j,iblock) = dlvb(i,j,iblock) +(SBCY(i,j,iblock)-BBCY(i,j,iblock))*ohbu(i,j,iblock)
!             ENDDO
!          ENDDO
!      ENDDO


! !$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
!     DO IBLOCK = 1, NBLOCKS_CLINIC
!       DO K = 1,KM
!          DO J = 1, JMT
!             DO I = 1, IMT
!                WKA (I,J,K,IBLOCK)= C0F * SQRT (UP (I,J,K,IBLOCK)* UP (I,J,K,IBLOCK) + VP (I, &
!                             J,K,IBLOCK)* VP (I,J,K,IBLOCK))
!             END DO
!          END DO
!       END DO
!    END DO
 
 
! !$OMP PARALLEL DO PRIVATE (IBLOCK,diff_u1,diff_v1,diff_u2,diff_v2)
!     DO IBLOCK = 1, NBLOCKS_CLINIC
!       DO K = 1,KM
!       DO J = 2, JMT-1
!          DO I = 2,IMT-1
! !lhl1204
! #ifdef CANUTO2010
! !       print*,akmU(i,j,k)
!             if (k==1) then
!                diff_v1 = SV (I,J,IBLOCK)* OD0*(1-AIDIF)
!                diff_u1 = SU (I,J,IBLOCK)* OD0*(1-AIDIF)
!             else
!                diff_v1= AKMU(I,J,K-1,IBLOCK)*(1-AIDIF)*(VP(I,J,K-1,IBLOCK)- VP(I,J,K,IBLOCK))* &
!                         ODZT(K)*VIV(I,J,K,IBLOCK)+ &
!                         (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK)*(1-AIDIF) &
!                       * (-SNLAT(I,J,IBLOCK)*UP(I,J,K-1,IBLOCK)*SAG+VP(I,J,K-1,IBLOCK)*CAG)
!                diff_u1= AKMU(I,J,K-1,IBLOCK)*(1-AIDIF)* (UP(I,J,K-1,IBLOCK)-UP(I,J,K,IBLOCK))*  &
!                         ODZT (K)*VIV (I,J,K,IBLOCK) + &
!                        (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK)*(1-AIDIF) &
!                       *(UP(I,J,K-1,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP(I,J,K-1,IBLOCK)* SAG)
!             end if
!             if (k==km) then
!                diff_v2= WKA (I,J,KM,IBLOCK)* ( - SNLAT (I,J,IBLOCK)* UP (I,J,KM,IBLOCK)        &
!                         * SAG + VP (I,J,KM,IBLOCK)* CAG)*(1-AIDIF)
!                diff_u2= WKA (I,J,KM,IBLOCK)* ( UP (I,J,KM,IBLOCK)* CAG + SNLAT (I,J,IBLOCK)    &
!                         * VP (I,J,KM,IBLOCK)* SAG)*(1-AIDIF)
!             else
!                diff_v2= AKMU(I,J,K,IBLOCK)*(1-AIDIF)*(VP(I,J,K,IBLOCK)- VP(I,J,K+1,IBLOCK))* &
!                         ODZT(K+1)*VIV(I,J,K+1,IBLOCK)+ &
!                         (1.0D0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK)*(1-AIDIF) &
!                       * (-SNLAT(I,J,IBLOCK)*UP(I,J,K,IBLOCK)*SAG+VP(I,J,K,IBLOCK)*CAG)
!                diff_u2= AKMU(I,J,K,IBLOCK)*(1-AIDIF)*(UP(I,J,K,IBLOCK)-UP(I,J,K+1,IBLOCK))*  & 
!                         ODZT(K+1)*VIV (I,J,K+1,IBLOCK) + &
!                        (1.0D0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK)*(1-AIDIF) &
!                       *(UP(I,J,K,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP(I,J,K,IBLOCK)* SAG)
!             end if
! #else
! !
!             call exit_licom(sigAbort,'The false mixing option')
! !
! #endif
! !lhl1204
!             DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + ODZP (K)* (diff_v1-diff_v2)
!             DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) + ODZP (K)* (diff_u1-diff_u2)
!         END DO
!       END DO
!       END DO
!    END DO
!
!    if (isc < 11 .and. mytid ==5) then
!       write(120+mytid,*) "DLU,  OK-3, SC=", isc, dlu(14,16,1,1),dlv(14,16,1,1)
!       write(120+mytid,*) "SU, SV,K-3, SC=", isc, su(14,16,1),sv(14,16,1)
!     end if
!     if (isc == 11 .and. mytid ==5) close(120+mytid)
!     deallocate(riu) 
 
!Yu   if (ierr /= 0) then
!Yu      write(6,*)'allocation error---tmp1,tmp2'
!Yu      stop
!Yu   end if
!     stop

      RETURN
      END SUBROUTINE READYC_debug