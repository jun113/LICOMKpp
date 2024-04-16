!     =================
      SUBROUTINE READYC
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
use POP_HaloMod
use POP_GridHorzMod
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
!use canuto_2010_mod
use canuto_mod
      IMPLICIT NONE
!      REAL(r8)  :: WKP (KMP1)
      INTEGER   :: IWK,n2, iblock,kmb,iwk1 !LPF20160729
      REAL(r8)  :: WK1 (KM-1) ,WK2 (KM-1), WK3 (KM-1),WK4(KM)
      REAL(r8)  :: WP1 (KM) ,WP2 (KM), WP3 (KM) !LPF20160728
      !REAL(r8)  :: WP1 (KM-1) ,WP2 (KM-1), WP3 (KM-1) !LPF20160728
      !REAL(r8)  :: WP1 (KM) ,WP2 (KM), WP3 (KM)
      REAL(r8)  :: WP4 (KM) ,WP5 (KM), WP6 (KM)
      !REAL(r8)  :: WP4 (KMM1) ,WP5 (KMM1), WP6 (KMM1)
      REAL(r8)  :: WP7 (KM) ,WP8(KM),tau_mag,mldtmp !LPF20160715
      !REAL(r8)  :: WP7 (KMM1) ,WP8(KM),tau_mag !LPF20160715
      REAL(r8)  :: WP9 ,WP10, WP11
      REAL(r8),dimension(IMT,JMT,KM,MAX_BLOCKS_CLINIC) :: WP12,WP13,riv1,riv2
      REAL(r8)  :: epsln,RKV,RKV1,ek0
      REAL(r8)  :: adv_x1,adv_x2,adv_x,adv_y1,adv_y2,adv_z,diff_u1,diff_u2,diff_v1,diff_v2
      REAL(r8)  :: dlux,dlvx,dluy,dlvy,dluz,dlvz,adv_z1,adv_z2,adv_z3,adv_z4
      real(r8)  :: hdvk(imt,jmt), hduk(imt,jmt), adv_uu(imt,jmt,km), adv_vv(imt,jmt,km)
      real(r8)  :: tq,sq
!LPF20160505
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
      integer ::  ErrorCode

       
#if (defined CANUTO)
      REAL(r8)  :: AIDIF
#endif

!YU 
      real (r8):: akt_back(km-1),aks_back(km-1),akm_back(km-1)
      !real (r8):: akt_back(2*km),aks_back(2*km),akm_back(2*km) !LPF20160728
      allocate(riu(imt,jmt,0:km,max_blocks_clinic),stat=ierr) 
!Yu   if (ierr /= 0) then
!Yu      write(6,*)'allocation error---riu'
!Yu      stop
!Yu   end if
 
      epsln = 1.0D-25
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JET
         DO I = 1,IMT
            H0BL (I,J,IBLOCK)= H0BF (I,J,IBLOCK)
            H0BF (I,J,IBLOCK)= H0 (I,J,IBLOCK)
         END DO
      END DO
   END DO
 
!lhl0711
#if (defined CANUTO)
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


!$OMP PARALLEL DO PRIVATE (IBLOCK,K)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         call ugrid_to_tgrid(wp12(:,:,k,iblock),up(:,:,k,iblock),iblock,k)
         call ugrid_to_tgrid(wp13(:,:,k,iblock),vp(:,:,k,iblock),iblock,k)
      END DO
   END DO
!
 
!$OMP PARALLEL DO PRIVATE (IBLOCK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KMM1
         DO J = 2, JMT
            DO I = 1,IMT-1
               riv1(i,j,k,iblock) = wp12 (I,J,K,iblock)*vit(i,j,k,iblock) - wp12 (I,J,K +1,iblock)*vit(i,j,k+1,iblock)
               riv2(i,j,k,iblock) = wp13 (I,J,K,iblock)*vit(i,j,k,iblock) - wp13 (I,J,K +1,iblock)*vit(i,j,k+1,iblock)
               s2t (i,j,k,iblock) =vit(i,j,k+1,iblock)*(riv1(i,j,k,iblock)*riv1(i,j,k,iblock)+riv2(i,j,k,iblock)*riv2(i,j,k,iblock))*ODZT(K+1)*ODZT(K+1)
#ifdef CANUTO  
               rit (i,j,k,iblock)= VIT (I,J,K +1,iblock)*rict(i,j,k,iblock)/(s2t(i,j,k,iblock)+epsln)
#else          
               rit (i,j,k,iblock)= rit (i,j,k,iblock,iblock) +VIT (I,J,K +1,iblock,iblock)* &
                                   rict(i,j,k,iblock,iblock)/(s2t(i,j,k,iblock,iblock)+epsln)
#endif
            END DO
         END DO
      END DO
   END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KMM1
         DO J = 1, JMT
            DO I = 1, IMT
               riv1(i,j,k,iblock) = UP (I,J,K,iblock) - UP (I,J,K +1,iblock)
               riv2(i,j,k,iblock) = VP (I,J,K,iblock) - VP (I,J,K +1,iblock)
               s2u (i,j,k,iblock) =viv(i,j,k+1,iblock)*(riv1(i,j,k,iblock)*riv1(i,j,k,iblock)+riv2(i,j,k,iblock)*riv2(i,j,k,iblock))*ODZT(K+1)*ODZT(K+1)
! calculate the shear square and the T component minus S component of Richardson Number
!lhl1204
               riu (i,j,k,iblock) = VIV (I,J,K +1,iblock)*ric (i,j,k,iblock)/(s2u(i,j,k,iblock)+epsln)
!lhl1204
            END DO
         END DO
      END DO
   END DO
!

!LPF20160505
#ifdef BCKMEX
!add the vertical diffusivity due to internal wave mixing based on Jochum(2009)
!$OMP PARALLEL
!$OMP WORKSHARE
             diff_back=0.0d0
             diff_back_sh=0.0d0
             diff_back_nh=0.0d0
!$OMP END WORKSHARE
!$OMP END PARALLEL
!$OMP PARALLEL DO PRIVATE(iblock,j,i)
   DO iblock = 1, nblocks_clinic
     DO J = 1, JMT
         DO I = 1, IMT
          if(tlat(i,j,iblock)<0.0) then
           diff_back_sh(i,j,iblock) =diff_back_coef_max*exp(-(0.4d0*(tlat(i,j,iblock)/DEGtoRAD+28.9))**2)
          else
           diff_back_nh(i,j,iblock) = diff_back_coef_max*exp(-(0.4d0*(tlat(i,j,iblock)/DEGtoRAD-28.9))**2)
          endif
          diff_back(i,j,iblock)=diff_back_eq+diff_back_sh(i,j,iblock)+diff_back_nh(i,j,iblock)
        if ( tlat(i,j,iblock) .lt. -10.0*degtorad ) then
          diff_back(i,j,iblock) = diff_back(i,j,iblock) + diff_back_coef
        elseif  ( abs(tlat(i,j,iblock)) .le. 10.0*degtorad ) then
          diff_back(i,j,iblock) = diff_back(i,j,iblock) + diff_back_coef*(abs(tlat(i,j,iblock)/DEGtoRAD)/10.0)**2
        else
          diff_back(i,j,iblock) = diff_back(i,j,iblock) + diff_back_coef
        endif
             !    if(mytid==0) write(*,*)'90-tlat(i,j,iblock)/DEGtoRAD,diff_back_sh=,diff_back_nh=',90-tlat(i,j,iblock)/DEGtoRAD,diff_back(i,j,iblock),&
             !       diff_back_sh(i,j,iblock),diff_back_nh(i,j,iblock) 
       ENDDO
     ENDDO
   ENDDO
!OMP END PARALLEL DO
#endif 
! BACKMX
!LPF20160505
       !if(mytid==0) write(*,*)',mytid,isc=',mytid,isc
#ifdef CANUTO
!
     amld = c0
     akmt = c0
     akmu = c0
     akm_back=c0
     akt_back=c0
     aks_back=c0


!
!$OMP PARALLEL DO PRIVATE (iblock,tq,sq,wp1,wp2,wp3,wp4,wp5,wp6,wp7,wp8,wp9,wp10,wp11,iwk,iwk1,wk1,wk2,wk3,wk4)
   do iblock = 1, nblocks_clinic
      DO J = 2, JMT-1
         DO I = 2, IMT-1

        if (VIT(I,J,1,iblock).gt.0.5) then
!
         wp1=0.D0
         wp2=0.D0
         wp3=0.D0
         wp4=0.D0
         wp5=0.D0
         wp6=0.D0
         wp7=0.D0
         wp8=0.D0
!
         do k=1,KMT(i,j,iblock) - 1
         !do k=1,KM - 1
            wp8(k)=-vit(i,j,k+1,iblock)*ZKP(k+1)*1.0d+2
            !wp8(k)=-vit(i,j,k+1,iblock)*ZKP(k+1)
         end do
!
         DO K = 1,KMT(i,j,iblock)-1
               !wp3(K)= VIT (I,J,K+1,iblock)* (pdensity (I,J,K,iblock)+(pdensity (I,J,K+1,iblock)-pdensity(I,J,K,iblock))* &
               !                DZP (K)/(DZP(K)+DZP(K+1)))*1.d-3 !lhl130125
               wp1(K)= VIT (I,J,K+1,iblock)* (AT (I,J,K,1,iblock)-(AT (I,J,K,1,iblock)-AT (I,J,K+1,1,iblock))* &
                               DZP (K)/(DZP(K)+DZP(K+1)))
               wp2(K)= VIT (I,J,K+1,iblock)* (AT (I,J,K,2,iblock)-(AT (I,J,K,2,iblock)-AT (I,J,K+1,2,iblock))* &
                               DZP (K)/(DZP(K)+DZP(K+1))) !for canuto2002
               !wp2(K)= VIT (I,J,K+1,iblock)* (AT (I,J,K,2,iblock)-(AT (I,J,K,2,iblock)-AT (I,J,K+1,2,iblock))* &
               !                DZP (K)/(DZP(K)+DZP(K+1)))*1.0D3+35.
               wp4(k)=vit(i,j,k+1,iblock)*RIT(i,j,k,iblock)
               wp5(k)=vit(i,j,k+1,iblock)*ricdtTmS(i,j,k,iblock) !for 
               !wp5(k)=vit(i,j,k+1,iblock)*RICDT(i,j,k,iblock)   !for canuto2010 !LPF20160728
               wp6(k)=vit(i,j,k+1,iblock)*S2T(i,j,k,iblock)
               wp7(k)=vit(i,j,k+1,iblock)*RICT(i,j,k,iblock)
               tq = wp1(k) - to(1)
               sq = (wp2(k)-35.0D0)*1.0D-3 - so(1)
               !wp3(k) = dens(tq,sq,1) !for canuto2010
               wp3(K)= VIT (I,J,K+1,iblock)*(pdensity (I,J,K,iblock)+&
                       (pdensity (I,J,K+1,iblock)-pdensity(I,J,K,iblock))* &
                               DZP (K)/(DZP(K)+DZP(K+1)))*1.d-3 !lhl130125 &LPF20160728
         END DO
         
        !!LPF20160728 
         wp9=vit(i,j,1,iblock)*USTAR(I,J,iblock)*1.0d+2
         wp10=vit(i,j,1,iblock)*BUOYTUR(I,J,iblock)*1.0d+4 
         wp11=vit(i,j,1,iblock)*BUOYSOL(I,J,iblock)*1.0d+4

#if ( defined CANUTOMIXOUT )
         wp10_canuto(I,J,iblock)=wp10     !yuzp-2016/12/8
         wp11_canuto(I,J,iblock)=wp11     !yuzp-2016/12/8
#endif
 
         IWK=KMT(I,J,iblock)-1 !need to be test
         IWK1=max(1,KMT(I,J,iblock)-1) !need to be test
        !!LPF20160728 
         tau_mag=USTAR(i,j,iblock)*USTAR(i,j,iblock)/OD0 !for canuto2010 !LPF20160715

!
!
!input ZKT in cm, AT(1) in C, AT(2) in (s-35)/1000., PDENSITY in g/cm^3
!        if (mytid ==0 .and. iday == 1 .and. isc == 0) then
!           write(180,*) i,j, iblock
!           write(180,*) wp1,wp2,wp3, wp8
!           write(180,*) wp4,wp5,wp6, wp7
!           write(180,*) wp9, wp10, wp11
!           write(180,*) 
!           write(180,*) (wp12(i,j,k,1),k=1,30)
!           write(180,*) (wp13(i,j,k,1),k=1,30)
!           write(180,*) 
!           write(180,*) (up(i,j,k,1),k=1,30)
!           write(180,*) (vp(i,j,k,1),k=1,30)
!           write(180,*) 
!           write(180,*) (up(i+1,j,k,1),k=1,30)
!           write(180,*) (vp(i+1,j,k,1),k=1,30)
!           write(180,*) 
!           write(180,*) (up(i,j-1,k,1),k=1,30)
!           write(180,*) (vp(i,j-1,k,1),k=1,30)
!           write(180,*) 
!           write(180,*) (up(i+1,j-1,k,1),k=1,30)
!           write(180,*) (vp(i+1,j-1,k,1),k=1,30)
!           write(180,*) 
!           write(180,*) (s2t(i,j,k,1),k=1,29)
!           close(180)
!        end if

#ifdef CANUTO2010
         if(mytid==0) write(*,*)"into canuto2010"
         call  canuto_2010_interface(wk1,    wk2,     wk3, wk4, amld(i,j,iblock) ,tau_mag,&
                                     wp1,    wp2,     wp3, wp4, wp5,              &
                                     wp7,    wp6,     ulat(i,j,iblock)/degtorad ,     wp8,& 
                                     kmt(i,j,iblock) ,i,j, iblock, isc)
#endif

!        if (mytid == 6 .and. isc ==41 .and. iday == 9 .and. i == 32 .and. j == 17) then
!          write(160,*) i,j,kmt(i,j,iblock)
!           write(160,*) wk1, wk2, wk3
!           write(160,*) wp1
!           write(160,*) wp2
!           write(160,*) wp3
!           write(160,*) wp4
!           write(160,*) wp5
!           write(160,*) wp6
!           write(160,*) wp7
!           write(160,*) wp8
!           close(160)
!        end if
!        if (isc ==41 .and. iday == 9 .and. i == 32 .and. j == 17) stop
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'testinput,iwk,iblock,isc',iwk,iblock,isc,vit(i,j,1,iblock) !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp8',iwk,wp8*100 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp1',wp1 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp2',wp2 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp3',wp3 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp4',wp4 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp5',wp5 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp6',wp6 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp7',wp7 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp9',wp9 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp10',wp10 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'wp11',wp11 !LPF20160728
        !if(mytid==32.and.i==21.and.j==22)write(160,*)'fcort',fcort(i,j,iblock) !LPF20160728
        CALL  TURB_2(wp8(1),wp1(1),wp2(1),wp3(1),&
!input RIT and RIDT no unit, S2T in 1/s^2, DFRICMX and DWNDMIX in cm^2/s
               wp4(1),wp5(1),wp6(1), DFRICMX*1.0d+4, DWNDMIX*1.0d+4,&
!output in cm^2/s, so 1d-4 should be multipled
              AKM_BACK(1),AKT_BACK(1),AKS_BACK(1),&
!input  RICT in 1/s^2 USTAR in cm/s, BUOYTUR,BUOYSOL in cm^2/s^3,FF in 1/s
              wp7(1),wp9,wp10,wp11,FCORT(i,j,iblock) ,& !OK
!output amld in cm, akmt,akh, and aks in cm^2/s
!               AMLD(I,J,iblock),AKMT(I,J,1,iblock),AKT(I,J,1,1,iblock),AKT(I,J,1,2,iblock),&
              mldtmp,WK1(1),WK2(1),WK3(1),&
!input int
              IWK,iwk1,KM-1,1,0,0,i,j,iblock) !OK!
              !IWK,kmt(I,J,iblock)-1,KM,1,0,0,i,j,iblock) !OK!
 !if(mytid==32.and.i==21.and.j==22)write(160,*)'wk1,isc',wk1,isc,imt,jmt !LPF20160728
 !if(mytid==32.and.i==21.and.j==22)write(160,*)'wk12',wk2 !LPF20160728
 !if(mytid==32.and.i==21.and.j==22)write(160,*)'wk13',wk3 !LPF20160728
 !if(mytid==32.and.i==21.and.j==22)write(160,*)'amld',amld(i,j,iblock) !LPF20160728

!change cm to meter in canuto2002
               AMLD(I,J,iblock)=mldtmp*1.0d-2 !LPF20160803 

#if ( defined TIDEMIX )
          
          AK_TIDE(I,J,:,IBLOCK)=0.0
          !LPF2016Nov12 IF (1./OHBT(I,J,IBLOCK).GT.SHELF_CUTOFF) THEN
          !IF (1./OHBT(I,J,IBLOCK).GT.SHELF_CUTOFF) THEN
          DO K = 1,INT(KMT(I,J,IBLOCK))-1
!for Indonesian Seas (-8.5~-1, 103~142)
!          IF (((tlon(i,j,iblock).gt.103.).and.(tlon(i,j,iblock).lt.142.))&
!          .and.((tlat(i,j,iblock).gt.-8.5).and.(tlat(i,j,iblock).lt.-1.))) THEN
!
!          WAVE_DIS(I,J,IBLOCK)=dmin1(WAVE_DIS(I,J,IBLOCK),max_wavedis)
!          AK_TIDE(I,J,K,iblock)=back_tidalmixing+(MIXING_EF*WAVE_DIS(I,J,IBLOCK)*FZ_TIDE(I,J,K,IBLOCK)&
!                  /DMAX1(RICDT(I,J,K,IBLOCK),1.D-8)/(WP3(K)+PO(K)+1000.0)*VIT(I,J,K+1,IBLOCK))
!          AK_TIDE(I,J,K,iblock)=DMIN1(AK_TIDE(I,J,K,iblock),max_tidalmixing)
!          AK_TIDE(I,J,K,iblock)=dmin1((mixing_ef*wave_dis(I,J,iblock)*fz_tide(I,J,K,iblock)&
!                  /dmax1(RICDT(I,J,int(kmt(i,j,iblock))-1,iblock),5.3265D-9)/&
!                  (wp3(int(kmt(i,j,iblock))-1)+po(int(kmt(i,j,iblock))-1)+1000.0)*vit(i,j,k+1,iblock))+1.D-5,5.D-3)
!          ELSE
!          AK_TIDE(I,J,K,iblock)=dmin1((local_mixing_fraction*mixing_ef*wave_dis(I,J,iblock)*fz_tide(I,J,K,iblock)&
!                  *RICDT(I,J,K,iblock)/(RICDT(I,J,K,iblock)+5.3275D-9)/RICDT(I,J,K,iblock)/&
!                  (wp3(K)+po(K)+1000.0)*vit(i,j,k+1,iblock))+1.D-5,5.D-3)
!             WAVE_DIS(I,J,IBLOCK)=dmin1(WAVE_DIS(I,J,IBLOCK),max_wavedis) !LPF2016Nov12
          !AK_TIDE(I,J,K,iblock)=back_tidalmixing+(LOCAL_MIXING_FRACTION*MIXING_EF*WAVE_DIS(I,J,IBLOCK)*FZ_TIDE(I,J,K,IBLOCK)&
          !        /DMAX1(RICDT(I,J,K,IBLOCK),1.D-8)/(WP3(K)+PO(K)+1000.0)*VIT(I,J,K+1,IBLOCK))
!          AK_TIDE(I,J,K,iblock)=back_tidalmixing+(LOCAL_MIXING_FRACTION*MIXING_EF*WAVE_DIS(I,J,IBLOCK)*FZ_TIDE(I,J,K,IBLOCK)&
!                  *RICDT(I,J,K,iblock)/(RICDT(I,J,K,IBLOCK)+5.3275D-9)/dmax1(RICDT(I,J,K,iblock),1.0d-8)/(WP3(K)*1000.0)*VIT(I,J,K+1,IBLOCK))

!          AK_TIDE(I,J,K,iblock)=back_tidalmixing+(LOCAL_MIXING_FRACTION*MIXING_EF*WAVE_DIS(I,J,IBLOCK)*FZ_TIDE(I,J,K,IBLOCK)&     !yuzp-2016/11/13
!                  /dmax1(ricdtTmS(I,J,K,iblock),1.0d-8)/(WP3(K)+PO(K)+1000.0)*VIT(I,J,K+1,IBLOCK))     !yuzp-2016/11/13

!          AK_TIDE1(I,J,K,iblock)=AK_TIDE(I,J,K,iblock)     !yuzp-2016/11/19


                  !*RICDT(I,J,K,iblock)/(RICDT(I,J,K,IBLOCK)+5.3275D-9)/dmax1(RICDT(I,J,K,iblock),1.0d-8)/(WP3(K)+PO(K)+1000.0)*VIT(I,J,K+1,IBLOCK))
!          AK_TIDE(I,J,K,iblock)=DMIN1(AK_TIDE(I,J,K,iblock),max_tidalmixing) !LPF2016Nov12
!                  /DMAX1(RICDT(I,J,INT(KMT(I,J,IBLOCK))-1,IBLOCK),1.D-8)/&
!                  /dmax1(RICDT(I,J,int(kmt(i,j,iblock))-1,iblock),5.3265D-9)/&
!                  *RICDT(I,J,int(kmt(i,j,iblock))-1,iblock)/(RICDT(I,J,int(kmt(i,j,iblock))-1,iblock)+5.3275D-9)/RICDT(I,J,int(kmt(i,j,iblock))-1,iblock)/&
!                 *fz_tide(I,J,K,iblock)/dmax1(RICDT(I,J,K,iblock),1.d-8)/(wp3(K)+po(K)+1000.0)*vit(i,j,k+1,iblock))+1.D-5,5.D-3)
!                 *fz_tide(I,J,K,iblock)/dmax1(RICDT(I,J,K,iblock),1.d-8)/wp3(K)*vit(i,j,k+1,iblock))+1.D-5,5.D-3)
!                 *fz_tide(I,J,K,iblock)/(dmax1(RICDT(I,J,K,iblock),0.0)+5.3275D-9)/wp3(K)*vit(i,j,k+1,iblock))+1.D-5,5.D-3)
!          ENDIF

          AK_TIDE(I,J,K,iblock)=back_tidalmixing+mixing_ef*local_mixing_fraction &
                                *WAVE_DIS(I,J,IBLOCK)*FZ_TIDE(I,J,K,IBLOCK) &
                               /(dmax1(rict(I,J,K,iblock),1.0d-8)*WP3(K)*1000.0)     !yuzp-2016/11/23--12/2

          AK_TIDE(I,J,K,iblock)=DMIN1(AK_TIDE(I,J,K,iblock),max_tidalmixing)   !LPF2016Nov12

          richardson(I,J,K,iblock)=rict(I,J,K,IBLOCK)     !yuzp-2016/11/13--12/2
          fztidal(I,J,K,iblock)=FZ_TIDE(I,J,K,IBLOCK)     !yuzp-2016/11/13
          wp3_tidal(I,J,K,iblock)=WP3(K)     !yuzp-2016/11/13

          END DO
          
          !fcor_canuto(I,J,iblock)=FCOR(I,J,IBLOCK)     !yuzp-2016/12/4
          !fcort_canuto(I,J,iblock)=FCORT(i,j,iblock)     !yuzp-2016/12/4

#if ( defined CANUTOMIXOUT )
          DO K = 1,KMT(i,j,iblock)-1
          wp1_canuto(I,J,K,iblock)=WP1(K)     !yuzp-2016/12/4
          wp2_canuto(I,J,K,iblock)=WP2(K)     !yuzp-2016/12/4
          wp3_canuto(I,J,K,iblock)=WP3(K)     !yuzp-2016/12/4
          wp4_canuto(I,J,K,iblock)=WP4(K)     !yuzp-2016/12/4
          wp5_canuto(I,J,K,iblock)=WP5(K)     !yuzp-2016/12/4
          wp6_canuto(I,J,K,iblock)=WP6(K)     !yuzp-2016/12/4
          wp7_canuto(I,J,K,iblock)=WP7(K)     !yuzp-2016/12/4
          wp8_canuto(I,J,K,iblock)=WP8(K)     !yuzp-2016/12/4

          wp12_canuto(I,J,K,iblock)=wp12(I,J,K,iblock)     !yuzp-2016/12/4
          wp13_canuto(I,J,K,iblock)=wp13(I,J,K,iblock)     !yuzp-2016/12/4

          wk4_canuto(I,J,K,iblock)=WK4(K)     !yuzp-2016/12/4

          END DO
#endif

!          K = 1,KMT(i,j,iblock)-1
!          alpha_canuto(I,J,K,iblock)=ALPHA(I,J,K,iblock)     !yuzp-2016/12/4
!          beta_canuto(I,J,K,iblock)=BETA(I,J,K,iblock)     !yuzp-2016/12/4
!          wp12_canuto(I,J,K,iblock)=wp12(I,J,K,iblock)     !yuzp-2016/12/4
!          wp13_canuto(I,J,K,iblock)=wp13(I,J,K,iblock)     !yuzp-2016/12/4
!          END DO

          DO K = INT(KMT(I,J,IBLOCK))-2,1,-1     !yuzp-2016/11/23
           AK_TIDE(I,J,K,iblock)=DMIN1(AK_TIDE(I,J,K,iblock),AK_TIDE(I,J,K+1,iblock))     !yuzp-2016/11/23
          END DO     !yuzp-2016/11/23

          !ENDIF !LPF2016Nov12
!      if (mytid.eq.master_task) then
!      write(6,*) kmt(i,j,1)
!      write(6,*) wave_dis(i,j,1)
!      write(6,*) RICDT(i,j,int(kmt(i,j,iblock))-1,1)
!      write(6,*) RICDT(i,j,int(kmt(i,j,iblock)),1)
!      write(6,*) dmax1(RICDT(i,j,:,1),1.d-8)
!      write(6,*) AK_TIDE(i,j,:,1)
!      write(6,*) wp3(:)+po(:)+1000.
!      write(6,*) wp3(:)
!      write(6,*) wp3(int(kmt(i,j,iblock))-1)+po(int(kmt(i,j,iblock))-1)+1000.
!      write(6,*) wp3(int(kmt(i,j,iblock)))+po(int(kmt(i,j,iblock)))+1000.
!      endif


#endif

         DO K = 1,KM-1
#if ( defined TIDEMIX )
!         AKMT(I,J,K,iblock)=WK1(K)+AK_TIDE(I,J,K,iblock)
         AKMT(I,J,K,iblock)=WK1(K)*1.0d-4+AK_TIDE(I,J,K,iblock)*5.0
         !AKMT(I,J,K,iblock)=WK1(K)*1.0d-4+AK_TIDE(I,J,K,iblock)*10.0
         AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+ (WK2(K)*1.0d-4+AK_TIDE(I,J,K,iblock))/float(NCC)
         AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+ (WK3(K)*1.0d-4+AK_TIDE(I,J,K,iblock))/float(NCC)

#if ( defined CANUTOMIXOUT )
         wk1_canuto(I,J,K,iblock)=WK1(K)     !yuzp-2016/12/4
         wk2_canuto(I,J,K,iblock)=WK2(K)     !yuzp-2016/12/4
         wk3_canuto(I,J,K,iblock)=WK3(K)     !yuzp-2016/12/4
#endif

!
#ifdef BCKMEX
          AKMT(I,J,K,iblock)=AKMT(I,J,K,iblock)+diff_back(i,j,iblock)*10.0*1.d-4 !10--Pr_number
          AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+diff_back(i,j,iblock)/float(NCC)*1.d-4
          AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+diff_back(i,j,iblock)/float(NCC)*1.d-4
#endif


!         IF (abs(ZKP(K+1)).GT.AMLD(I,J,iblock)) THEN
!         AKMT(I,J,K,iblock)=dmin1(AKMT(I,J,K,iblock),1.5D-1)
!         AKMT(I,J,K,iblock)=dmax1(AKMT(I,J,K,iblock),0.1D-3)
!         AKT(I,J,K,1,iblock)=dmin1(AKT(I,J,K,1,iblock),2.4D-1)
!         AKT(I,J,K,1,iblock)=dmax1(AKT(I,J,K,1,iblock),0.2D-4)
!         !AKT(I,J,K,1,iblock)=dmax1(AKT(I,J,K,1,iblock),0.5D-4)
!         AKT(I,J,K,2,iblock)=dmin1(AKT(I,J,K,2,iblock),2.4D-1)
!         AKT(I,J,K,2,iblock)=dmax1(AKT(I,J,K,2,iblock),0.2D-4)
!         !AKT(I,J,K,2,iblock)=dmax1(AKT(I,J,K,2,iblock),0.5D-4)
!         ELSE
!         AKMT(I,J,K,iblock)=dmax1(AKMT(I,J,K,iblock),1.D-2)
!         AKT(I,J,K,1,iblock)=dmax1(AKT(I,J,K,1,iblock),1.D-3)
!         AKT(I,J,K,2,iblock)=dmax1(AKT(I,J,K,2,iblock),1.D-3)
!         ENDIF
#else
         AKMT(I,J,K,iblock)=WK1(K)*1.0d-4
         AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+ WK2(K)*1.0d-4/float(NCC)
         AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+ WK3(K)*1.0d-4/float(NCC)

#ifdef BCKMEX
          AKMT(I,J,K,iblock)=AKMT(I,J,K,iblock)+diff_back(i,j,iblock)*10.0*1.d-4 !10--Pr_number
          AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+diff_back(i,j,iblock)/float(NCC)*1.d-4
          AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+diff_back(i,j,iblock)/float(NCC)*1.d-4
#endif
!         IF (abs(ZKP(K+1)).GT.AMLD(I,J,iblock)) THEN
!         AKMT(I,J,K,iblock)=dmin1(AKMT(I,J,K,iblock),1.5D-1)
!         AKMT(I,J,K,iblock)=dmax1(AKMT(I,J,K,iblock),0.1D-3)
!         AKT(I,J,K,1,iblock)=dmin1(AKT(I,J,K,1,iblock),2.4D-1)
!         AKT(I,J,K,1,iblock)=dmax1(AKT(I,J,K,1,iblock),0.2D-4)
!         !AKT(I,J,K,1,iblock)=dmax1(AKT(I,J,K,1,iblock),0.5D-4)
!         AKT(I,J,K,2,iblock)=dmin1(AKT(I,J,K,2,iblock),2.4D-1)
!         AKT(I,J,K,2,iblock)=dmax1(AKT(I,J,K,2,iblock),0.2D-4)
!         !AKT(I,J,K,2,iblock)=dmax1(AKT(I,J,K,2,iblock),0.5D-4)
!         ELSE
!         AKMT(I,J,K,iblock)=dmax1(AKMT(I,J,K,iblock),1.D-2)
!         AKT(I,J,K,1,iblock)=dmax1(AKT(I,J,K,1,iblock),1.D-3)
!         AKT(I,J,K,2,iblock)=dmax1(AKT(I,J,K,2,iblock),1.D-3)
!         ENDIF
#endif
         END DO


!               k=6
!              if (mytid == 7  .and. j == 25 .and. i == 18 ) then
!                 write(229,*) i,j,k
!                 write(229,*) wp8(k), wp1(k), wp2(k), wp3(k), wp4(k), wp5(k), wp6(k), DFRICMX, DWNDMIX
!                 write(229,*) akm_back(k), akt_back(k), aks_back(k)
!                 write(229,*) wp7(k), wp9, wp10, wp11, FCORT(i,j,1)
!                 write(229,*) amld(i,j,1), wk1(k), wk2(k), wk3(k), iwk, kmt(i,j,1)-1, km
!                 write(229,*) akmt(i,j,k,1)
!              end if
!
        endif
         END DO
      END DO
   END DO
!  stop
!!
! calculate the vertical mixing on U-grid

!$OMP PARALLEL DO PRIVATE (IBLOCK)
   do iblock = 1, nblocks_clinic
      DO K = 1,KMM1
         call tgrid_to_ugrid(akmu(:,:,k,iblock), akmt(:,:,k,iblock),iblock)
         do j= 1,jmt
         do i= 1,imt
            akmu(i,j,k,iblock) = akmu(i,j,k,iblock)* viv(i,j,k+1,iblock)
         end do
         end do
      END DO
   end do
!lhl241204
#endif


         if (.not. allocated(rict_ref)) allocate(rict_ref(imt,jmt,max_blocks_clinic))

!LPF20160515
   call POP_HaloUpdate(amld , POP_haloClinic, POP_gridHorzLocCenter,&
                         POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!LPF20160515

        !   rict_ref = rict(:,:,14,:) !why 14layer
   DO IBLOCK = 1,NBLOCKS_CLINIC
       DO J = 1,JMT
          DO I = 1, IMT
             rict_ref(i,j,iblock) = rict(i,j,14,iblock) !why 14layer
            DO K = 1,KMM1
            IF (abs(ZKP(K+1)).GT.AMLD(I,J,iblock)) THEN
            rict_ref(I,J,IBLOCK) = rict(I,J,K,IBLOCK)
               exit
            ENDIF
            END DO
         END DO
      END DO
   END DO



!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM: ZONAL COMPONENT
!---------------------------------------------------------------------
      CALL UPWELL (U,V,H0)
      
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
     dlu = c0
     dlv = c0
     wka = c0

!$OMP PARALLEL DO PRIVATE (IBLOCK,K)
   do iblock = 1, nblocks_clinic
       DO K = 1,KM
          call tgrid_to_ugrid(wka(:,:,k,iblock),ws(:,:,k,iblock), iblock)
       END DO
   end do

!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERMS
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (IBLOCK)
    DO IBLOCK = 1, NBLOCKS_CLINIC
        call advection_momentum(u(:,:,:,iblock),v(:,:,:,iblock),wka(:,:,:,iblock), &
                                dlu(:,:,:,iblock),dlv(:,:,:,iblock),iblock)
!       dlu(:,:,:,IBLOCK)=adv_uu
!       dlv(:,:,:,IBLOCK)=adv_vv
   END DO
!
!---------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL VISCOSITY
!---------------------------------------------------------------------
 
#if ( defined SMAG)
!
         CALL SMAG2 (K)
!
#if (defined SMAG_FZ )
#else
#endif
#else
#if (defined BIHAR)

!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,HDUK,HDVK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
          call hdiffu_del4(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
          do j = 3, jmt-2
          do i = 3, imt-2
             dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
             dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
          end do
          end do
      END DO
   END DO

#else
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,hduk,hdvk,K)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call hdiffu_del2(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
         DO J = 3, JMT-2
            DO I = 3,IMT-2
               dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
               dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
            END DO
         END DO
      END DO
   END DO
 
#endif
#endif
 
      allocate(dlub(imt,jmt,max_blocks_clinic),dlvb(imt,jmt,max_blocks_clinic),stat=ierr)
!---------------------------------------------------------------------
!     VERTICAL INTEGRATION
!---------------------------------------------------------------------
      CALL VINTEG (DLU,DLUB)
      CALL VINTEG (DLV,DLVB)
!
!$OMP PARALLEL DO PRIVATE (iblock,J,I,kmb) 
      do iblock = 1, nblocks_clinic
         DO J = 2, jmt-1
            DO I = 2,imt-1
               kmb = kmu(i,j,iblock)
              if (kmb>=1) then !LPF20170621
               SBCX(I,J,IBLOCK) = SU (I,J,IBLOCK)* OD0
               SBCY(I,J,IBLOCK) = SV (I,J,IBLOCK)* OD0
               BBCX(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+  &
                                          VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))   &
                          *(UP(I,J,kmb,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP (I,J,kmb,IBLOCK)*SAG)
               BBCY(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+  &  
                                          VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))&
                          *(-SNLAT(I,J,IBLOCK)*UP(I,J,kmb,IBLOCK)*SAG+VP(I,J,kmb,IBLOCK)*CAG)
              else
                SBCX(I,J,IBLOCK)=0.0D0
                SBCY(I,J,IBLOCK)=0.0D0
                BBCX(I,J,IBLOCK)=0.0D0
                BBCY(I,J,IBLOCK)=0.0D0
              endif
              dlub(i,j,iblock) = dlub(i,j,iblock) +(SBCX(i,j,iblock)-BBCX(i,j,iblock))*ohbu(i,j,iblock)
              dlvb(i,j,iblock) = dlvb(i,j,iblock) +(SBCY(i,j,iblock)-BBCY(i,j,iblock))*ohbu(i,j,iblock)
            ENDDO
         ENDDO
     ENDDO


!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1, IMT
               WKA (I,J,K,IBLOCK)= C0F * SQRT (UP (I,J,K,IBLOCK)* UP (I,J,K,IBLOCK) + VP (I, &
                            J,K,IBLOCK)* VP (I,J,K,IBLOCK))
            END DO
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,diff_u1,diff_v1,diff_u2,diff_v2)
    DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
      DO J = 2, JMT-1
         DO I = 2,IMT-1
!lhl1204
#ifdef CANUTO
!       print*,akmU(i,j,k)
            if (k==1) then
               diff_v1 = SV (I,J,IBLOCK)* OD0*(1-AIDIF)
               diff_u1 = SU (I,J,IBLOCK)* OD0*(1-AIDIF)
            else
               diff_v1= AKMU(I,J,K-1,IBLOCK)*(1-AIDIF)*(VP(I,J,K-1,IBLOCK)- VP(I,J,K,IBLOCK))* &
                        ODZT(K)*VIV(I,J,K,IBLOCK)+ &
                        (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK)*(1-AIDIF) &
                      * (-SNLAT(I,J,IBLOCK)*UP(I,J,K-1,IBLOCK)*SAG+VP(I,J,K-1,IBLOCK)*CAG)
               diff_u1= AKMU(I,J,K-1,IBLOCK)*(1-AIDIF)* (UP(I,J,K-1,IBLOCK)-UP(I,J,K,IBLOCK))*  &
                        ODZT (K)*VIV (I,J,K,IBLOCK) + &
                       (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK)*(1-AIDIF) &
                      *(UP(I,J,K-1,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP(I,J,K-1,IBLOCK)* SAG)
            end if
            if (k==km) then
               diff_v2= WKA (I,J,KM,IBLOCK)* ( - SNLAT (I,J,IBLOCK)* UP (I,J,KM,IBLOCK)        &
                        * SAG + VP (I,J,KM,IBLOCK)* CAG)*(1-AIDIF)
               diff_u2= WKA (I,J,KM,IBLOCK)* ( UP (I,J,KM,IBLOCK)* CAG + SNLAT (I,J,IBLOCK)    &
                        * VP (I,J,KM,IBLOCK)* SAG)*(1-AIDIF)
            else
               diff_v2= AKMU(I,J,K,IBLOCK)*(1-AIDIF)*(VP(I,J,K,IBLOCK)- VP(I,J,K+1,IBLOCK))* &
                        ODZT(K+1)*VIV(I,J,K+1,IBLOCK)+ &
                        (1.0D0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK)*(1-AIDIF) &
                      * (-SNLAT(I,J,IBLOCK)*UP(I,J,K,IBLOCK)*SAG+VP(I,J,K,IBLOCK)*CAG)
               diff_u2= AKMU(I,J,K,IBLOCK)*(1-AIDIF)*(UP(I,J,K,IBLOCK)-UP(I,J,K+1,IBLOCK))*  & 
                        ODZT(K+1)*VIV (I,J,K+1,IBLOCK) + &
                       (1.0D0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK)*(1-AIDIF) &
                      *(UP(I,J,K,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP(I,J,K,IBLOCK)* SAG)
            end if
#else
!
            call exit_licom(sigAbort,'The false mixing option')
!
#endif
!lhl1204
            DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + ODZP (K)* (diff_v1-diff_v2)
            DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) + ODZP (K)* (diff_u1-diff_u2)
        END DO
      END DO
      END DO
   END DO
!
      RETURN
      END SUBROUTINE READYC
 
 
