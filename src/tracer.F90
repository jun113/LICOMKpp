!  CVS: $Id: tracer.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE TRACER
!     =================
#include <def-undef.h>
use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use work_mod
use dyn_mod
use isopyc_mod
use forc_mod
use pmix_mod
use msg_mod
use smuvh
use advection
use blocks
use domain
use LICOM_Error_mod
use gather_scatter
use distribution
use buf_mod
#ifdef BIHAR
use hmix_del4
#else
use hmix_del2
#endif
      IMPLICIT NONE
 
      integer     :: n2, iblock,kt, irec
      REAL(r8)    :: AIDIF,C2DTTS,AA,FAW,FIW,ALF,RNCC,ABC,fil_lat1,fil_lat2
      REAL(r8)    :: HDTK(imt,jmt), adv_tt(imt,jmt,km),ek0 ,ek1
      REAL (r8)   :: DT2K(imt,jmt)
#ifdef SSSNORM
!lhl20130913
      REAL(r8)    :: ERR_norm1,ERR_norm2
!lhl20130913
      real (r8),dimension(:,:,:),allocatable::temp11
      real (r8),dimension(:,:),allocatable::temp12 
#endif

!Xiao Chan (Hereinafter XC for short)
      real(r8)    :: LAMDA(imt,jmt,km,max_blocks_clinic),wt1,wt2,adv_y,adv_x,adv_z,adv_x1,adv_x2
      real(r8)    :: LAMDA1(km)
      type (block) :: this_block          ! block information for current block
#ifdef LPFDIAG
      REAL(r8)    :: outsum 
#endif   
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::VTL_ori !for output dt diffusion
!XC

#ifdef LICOM_ENABLE_TEST_TIME
#undef  LICOM_ENABLE_TEST_TRACER
#define LICOM_ENABLE_TEST_TRACER
#endif
 
!---------------------------------------------------------------------
!     SET LOCAL CONSTANT
!---------------------------------------------------------------------
      ! deallocate(dlu,dlv,gg,ric,rict)
      if (allocated(dlu)) then
         deallocate(dlu)
      end if
      if (allocated(dlv)) then
         deallocate(dlv)
      end if
      if (allocated(gg)) then
         deallocate(gg)
      end if
      if (allocated(ric)) then
         deallocate(ric)
      end if
      if (allocated(rict)) then
         deallocate(rict)
      end if
 
      ! allocate(stf(imt,jmt,max_blocks_clinic),tf(imt,jmt,km,max_blocks_clinic))
      ! allocate(wkb(imt,jmt,km,max_blocks_clinic),wkc(imt,jmt,km,max_blocks_clinic),wkd(imt,jmt,km,max_blocks_clinic))
      !wjl 20231123
      if (.not. allocated(stf)) then
         allocate(stf(imt,jmt,max_blocks_clinic))
      end if
      if (.not. allocated(tf)) then
         allocate(tf(imt,jmt,km,max_blocks_clinic))
      end if
      if (.not. allocated(wkb)) then
         allocate(wkb(imt,jmt,km, max_blocks_clinic))
      end if
      if (.not. allocated(wkc)) then
         allocate(wkc(imt,jmt,km, max_blocks_clinic))
      end if
      if (.not. allocated(wkd)) then
         allocate(wkd(imt,jmt,km, max_blocks_clinic))
      end if
!lhl      AIDIF = 0.5*FLOAT(ISOP)
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
       fil_lat1=63.0D0
       fil_lat2=63.0D0
!      GAMMA = 0.0D0

      AIDIF = 0.5D0
 
      LAMDA=1.0D0/(15.D0*86400.D0)
     do k=1,45
      LAMDA1(k)=1.0D0/(5.D0*86400.D0)
     enddo
     do k=46,km
      LAMDA1(k)=1.0D0/(5.D0*86400.D0)/(k*k/5.0/5.0)
     enddo

      RNCC = 1.0D0/ FLOAT (NCC)
!
!
      if ( adv_tracer(1:8) == 'centered' )  then
      IF (IST >= 1)THEN
         C2DTTS = DTS *2.0D0
         AA = 0.5D0
      ELSE
         C2DTTS = DTS
         AA = 0.0D0
      END IF
      else if ( adv_tracer(1:5) == 'tspas') then
         C2DTTS = DTS
         AA = 0.5D0
      else 
        call exit_licom(sigAbort,'The false advection option for tracer')
      end if
 
!  DO iblock = 1, nblocks_clinic
!     DO K = 1, 23  
!     DO J = 1, JMT
!        DO I = 1,IMT
!           if ( tlat(i,j,iblock) > 45.0D0*degtorad .and. tlat(i,j,iblock) < &
!                75.D0*degtorad .and. tlon(i,j,iblock) > -100.0D0*degtorad .and.  &
!                tlon(i,j,iblock) <  10.0D0*degtorad) then
!              AKT(I,J,K,1,iblock)= AKT (I,J,K,1,iblock) +8.0D-3
!              AKT(I,J,K,2,iblock)= AKT (I,J,K,2,iblock) +8.0D-3
!           end if
!        END DO
!     END DO
!     END DO
!  END DO

!---------------------------------------------------------------------
!     PREPARATION FOR VERTICAL ADVECTIVE TERM
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO iblock = 1, nblocks_clinic
      DO J = 1, JMT
         DO I = 1,IMT
            H0F (I,J,iblock)= H0F (I,J,iblock)* ONBC
         END DO
      END DO
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO iblock = 1, nblocks_clinic
      DO J = 1, JMT
         DO I = 1,IMT
            STF (I,J,IBLOCK)= AA * H0F (I,J,IBLOCK) + (1.0D0- AA)* H0L (I,J,IBLOCK)
         END DO
      END DO
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               UTF (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK)* ONCC
               VTF (I,J,K,IBLOCK)= VTF (I,J,K,IBLOCK)* ONCC
            END DO
         END DO
      END DO
  END DO

 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               WKD (I,J,K,IBLOCK)= AA * UTF (I,J,K,IBLOCK) + (1.0D0- AA)* UTL (I,J,K,IBLOCK)
               WKB (I,J,K,IBLOCK)= AA * VTF (I,J,K,IBLOCK) + (1.0D0- AA)* VTL (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO

      CALL UPWELL (WKD,WKB,STF)
 
#if (defined NODIAG)
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
 
!-----------------------------------------------------------------------
!     PREPARATION FOR ISOPYCNAL DIFFUSION & ADVECTION
!-----------------------------------------------------------------------
#if (defined ISO)
      CALL ISOPYC
#endif
 
!@@@  COMPUTING DIFFUSION COEFFICIENT
    
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 2, JMT-1
            DO I = 2, IMT-1
#if (defined ISO)
               WKC (I,J,K,IBLOCK) = AHV + AHISOP(i,j,iblock) * K3 (I,K,J,3,IBLOCK)
#else
               WKC (I,J,K,IBLOCK) = AHV
#endif
!
#if (!defined CANUTO)
            AKT(I,J,K,1,IBLOCK)=WKC(I,J,K,IBLOCK)
            AKT(I,J,K,2,IBLOCK)=WKC(I,J,K,IBLOCK)
#endif
!
            END DO
         END DO
      END DO
   END DO
 
!
!-----------------------------------------------------------------------
!     SOLVE FOR ONE TRACER AT A TIME
!-----------------------------------------------------------------------
!     NTRA = 1 => TEMPERATURE
!     NTRA = 2 => SALINITY
 
      DO N = 1,NTRA
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM 
!---------------------------------------------------------------------

tf=0.0D0 !LPF20160829


!$OMP PARALLEL DO PRIVATE (IBLOCK,adv_tt)
   do iblock = 1, nblocks_clinic
      adv_tt = 0.0_r8
      this_block = get_block(blocks_clinic(iblock),iblock)
      call advection_tracer(wkd(:,:,:,iblock),wkb(:,:,:,iblock),ws(:,:,:,iblock), &
                            at(:,:,:,n,iblock),adv_tt,iblock,N,this_block)
                           !ax(:,:,:,N,iblock),ay(:,:,:,N,iblock),az(:,:,:,N,iblock)) !LPF20160811
      do k=1, km
      do j =3, jmt-2
      do i =3, imt-2
         tf(i,j,k,iblock) = adv_tt(i,j,k)*vit(i,j,k,iblock) 
!        tf(i,j,k,iblock) = adv_tt(i,j,k)*vit(i,j,k,iblock) +  &
!                          0.1*gamma*(restore(i,j,k,n,iblock)- atb(i,j,k,n,iblock))
!ZWP20170210 value before gamma oraiginal = 0.1
      end do
      end do
      end do
   end do
    
#ifdef CANUTO      
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 2, JMT-1
            DO I = 2, IMT-1
                 IF (AKT(I,J,1,N,IBLOCK).lt.dwndmix) AKT(I,J,1,N,IBLOCK)=dwndmix !LPF20160409 
#if (defined ISO)       
!#if ( defined TIDEMIX )
!               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK)+AK_TIDE(I,J,K,IBLOCK) + AHISOP(i,j,iblock) * K3 (I,K,J,3,IBLOCK) 
!#else
               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK) + AHISOP(i,j,iblock) * K3 (I,K,J,3,IBLOCK) 
!#endif
#else            
!#if ( defined TIDEMIX )
!               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK)+AK_TIDE(I,J,K,IBLOCK) 
!#else
               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK)
!#endif
#endif         
            END DO
         END DO
      END DO
   END DO
!
#endif

!-----------------------------------------------------------------------
!     COMPUTE THE ISOPYCNAL/DIPYCNAL MIXING
!-----------------------------------------------------------------------
!     XZ AND YZ ISOPYCNAL DIFFUSIVE FLUX ARE SOLVED EXPLICITLY;
!     WHILE ZZ COMPONENT WILL BE SOLVED IMPLICITLY.
 
#ifdef LPFDIAG
!LPF20160514
          if (mytid ==0) WRITE (16,*)'before ISOFLUX'
          call outsumglobal(TF,outsum)
          if (mytid ==0) WRITE (16,*)'outsum',outsum 
!LPF20160514
#endif
 
#if (defined ISO)
         CALL ISOFLUX (N)

#ifdef LPFDIAG
!LPF20160514
          if (mytid ==0) WRITE (16,*)'after ISOFLUX'
          call outsumglobal(TF,outsum)
          if (mytid ==0) WRITE (16,*)'outsum',outsum
!LPF20160514
#endif

#else 
 
#if ( defined SMAG)
         CALL SMAG3
#else
!-----------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL DIFFUSION TERMS:  ZONAL and Meridional COMPONENTS
!-----------------------------------------------------------------------
 
#if (defined BIHAR)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,hdtk)
   do iblock = 1, nblocks_clinic
   this_block = get_block(blocks_clinic(iblock),iblock)
   do k =1, km
!20200221 LPF-yyq
       call hdifft_del4(k,DT2K,HDTK,ATB(:,:,k,n,iblock),this_block)
       !call hdifft_del4(k,HDTK,ATB(:,:,k,n,iblock),this_block)
!20200221 LPF-yyq
      do j= 3, jmt-2
      do i= 3, imt-2
           TF (I,J,K,IBLOCK) = TF (I,J,K,IBLOCK) + HDTK(I,J)
           dx(i,j,k,N,iblock)=HDTK(i,j) !LPF20160822
      end do
      end do
           !dx(:,:,k,N,iblock)=HDTK(:,:) !LPF20160822
   end do
   end do
!
#else
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,hdtk)
   do iblock = 1, nblocks_clinic
   this_block = get_block(blocks_clinic(iblock),iblock)
   do k =1, km
      call hdifft_del2(k,HDTK,ATB(:,:,k,n,iblock),this_block)
      do j= 3, jmt-2
      do i= 3, imt-2
           TF (I,J,K,IBLOCK) = TF (I,J,K,IBLOCK) + HDTK(I,J)
           dx(i,j,k,N,iblock)=HDTK(i,j) !LPF20160822
      end do
      end do
           !dx(:,:,k,N,iblock)=HDTK(:,:) !LPF20160822
   end do
   end do
!
#endif
#endif
#endif

#ifdef LPFDIAG
!LPF20160514
          if (mytid ==0) WRITE (16,*)'after1 ISOFLUX'
          call outsumglobal(TF,outsum)
          if (mytid ==0) WRITE (16,*)'outsum',outsum 
!LPF20160514
#endif

!-----------------------------------------------------------------------
!     VERTICAL COMPONENT
!-----------------------------------------------------------------------

         IF (N == 1)THEN
#if (defined SOLAR)
!     SOLAR SHORTWAVE PENETRATION
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
    DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K=2,KM-1
        DO J=3, JMT-2
        DO I=3,IMT-2
            wt1= SWV (I,J,IBLOCK)*pen(k-1)*VIT(I,J,K,IBLOCK)
            wt2= SWV (I,J,IBLOCK)*pen(k)*VIT(I,J,K+1,IBLOCK)
            TF (I,J,K,IBLOCK)= TF(I,J,K,IBLOCK)+(wt1-wt2)*ODZP(K)
!
!lhl0105           penetrate(i,j,k)=(wt1-wt2)*ODZP(k)
           penetrate(i,j,k,iblock)=(wt1-wt2)*ODZP(k)
!            penetrate(i,j,k,IBLOCK)= wt2*ODZP(k)
!
        END DO
        END DO
        END DO
     END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
     DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J=3, JMT-2
        DO I=3,IMT-2
           wt1= SWV (I,J,IBLOCK)*pen(1)*VIT(I,J,2,IBLOCK)
           wt2= SWV (I,J,IBLOCK)*pen(km-1)*VIT(I,J,km,IBLOCK)
           TF(I,J,1,IBLOCK)=TF(I,J,1,IBLOCK)-ODZP(1)*wt1
           TF(I,J,km,IBLOCK)=TF(I,J,km,IBLOCK)+ODZP(km)*wt2
!
            penetrate(i,j, 1,IBLOCK)= -wt1*ODZP( 1) !LPF20160824
            penetrate(i,j,km,IBLOCK)= wt2*ODZP(km) !LPF20160828
            !penetrate(i,j,km,IBLOCK)= -wt2*ODZP(km) !LPF20160824
!
        END DO
        END DO
     END DO
#endif


!==================================
!From here ie. the code about SOLARCHLORO added by linpf
!==================================

#if (defined SOLARCHLORO)
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K=2,KM-1
        DO J=3,JMT-2
        DO I=3,IMT-2
            wt1= SWV(I,J,IBLOCK)*pen_chl(I,J,K-1,IBLOCK)*VIT(I,J,K,IBLOCK)
            wt2= SWV(I,J,IBLOCK)*pen_chl(I,J,K,IBLOCK)*VIT(I,J,K+1,IBLOCK)
            TF (I,J,K,IBLOCK)= TF(I,J,K,IBLOCK)+(wt1-wt2)*ODZP(K)
!
            penetrate(I,J,K,IBLOCK)=(wt1-wt2)*ODZP(K)
!
        END DO
        END DO
      END DO
  END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J=3,JMT-2
        DO I=3,IMT-2
           wt1= SWV(I,J,iblock)*pen_chl(I,J,1,iblock)*VIT(I,J,2,iblock)
           wt2= SWV(I,J,iblock)*pen_chl(I,J,km-1,iblock)*VIT(I,J,km,iblock)
           TF(I,J,1,iblock)=TF(I,J,1,iblock)-ODZP(1)*wt1
           TF(I,J,km,iblock)=TF(I,J,km,iblock)+ODZP(km)*wt2
!
            penetrate(I,J, 1,iblock)=-wt1*ODZP(1)
            penetrate(I,J,km,iblock)= wt2*ODZP(km)
!
        END DO
        END DO
  END DO
#endif

!==================================
!above code about SOLARCHLORO added by linpf
!==================================
         END IF
 
!     EDDY-DIFFUSION
 
 
        wt1=0
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K=2,KM-1
           DO J=3,JMT-2
           DO I=3,IMT-2
              wt1= WKC(I,J,K-1,IBLOCK)*(ATB(I,J,K-1,N,IBLOCK)-ATB(I,J,K,N,IBLOCK))*ODZT(K)*VIT(I,J,K,IBLOCK)
              wt2= WKC(I,J,K,IBLOCK)*(ATB(I,J,K,N,IBLOCK)-ATB(I,J,K+1,N,IBLOCK))*ODZT(K+1)*VIT(I,J,K+1,IBLOCK)
              TF (I,J,K,IBLOCK)= TF(I,J,K,IBLOCK)+ODZP(K)*(wt1-wt2)*(1.0D0-AIDIF)
!
              dz(i,j,k,N,IBLOCK)=ODZP(K)*(wt1-wt2)*(1.0D0-AIDIF)
!
           END DO
           END DO
        END DO
  END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
       DO J = 3,JMT-2
       DO I = 3,IMT-2
           wt1= WKC(I,J,1,IBLOCK)*(ATB(I,J,1,N,IBLOCK)-ATB(I,J,2,N,IBLOCK))*ODZT(2)*VIT(I,J,2,IBLOCK)
           wt2= WKC(I,J,km-1,IBLOCK)*(ATB(I,J,km-1,N,IBLOCK)-ATB(I,J,km,N,IBLOCK))*ODZT(km)*VIT(I,J,km,IBLOCK)
           TF(I,J,1,IBLOCK)=TF(I,J,1,IBLOCK)-ODZP(1)*wt1*(1.0D0-AIDIF)
           TF(I,J,km,IBLOCK)=TF(I,J,km,IBLOCK)+ODZP(km)*wt2*(1.0D0-AIDIF)
!
              dz(i,j, 1,N,IBLOCK)=-ODZP( 1)*wt1*(1.0-AIDIF)
              dz(i,j,km,N,IBLOCK)= ODZP(km)*wt2*(1.0-AIDIF)
!
       END DO
       END DO
  END DO        
  
          !if(mytid==100) write(16,*) 'TF,0-0',TF(21,10,1,1)
 !needto update the halo !LPF20150914 

!-----------------------------------------------------------------------
!     SET NEWTONIAN SURFACE BOUNDARY CONDITION
!-----------------------------------------------------------------------
 
         IF (N == 2)THEN
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)    
        DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = 1, JMT
               DO I = 1, IMT
                  IF (KMT(I,J,IBLOCK) > 0)THEN
#ifdef COUP
                    if (boundary_restore == 2) then
!lhl20130914         
! under seaice 50m/30days, all 50m/4yr (GAMMA=1/30days)
                    STF (I,J,IBLOCK) = SSF(I,J,IBLOCK)/ODZP(1)&
                      +(GAMMA*50.*ODZP(1))*(SSS(I,J,IBLOCK) &
                      -ATB (I,J,1,2,IBLOCK))*ifrac(i,j,IBLOCK)/ODZP(1)&
                      +(GAMMA*30./365./4.*50.*ODZP(1))*(SSS(I,J,IBLOCK)-ATB (I,J,1,2,IBLOCK))/ODZP(1)
!lhl20130914
!            STF (I,J,IBLOCK) = SSF(I,J,IBLOCK)/ODZP(1) + GAMMA * (RESTORE (I,J,1,N,iblock) - ATB (I,J,1,2,IBLOCK))/ODZP(1)
                    else !for full couple run
                         STF (I,J,IBLOCK) = SSF(I,J,IBLOCK)/ODZP(1)
                    end if
#else
#ifdef FRC_CORE
                     STF (I,J,IBLOCK) = (fresh(i,j,IBLOCK)*34.7*OD0*1.0D-3 & 
                     !STF (I,J,IBLOCK) = (fresh(i,j,IBLOCK)*34.7/1000.0/1000.0& 
                     !STF (I,J,IBLOCK) = (fresh(i,j)*35.0/1000.0/1000.0& !LPF20160817
                   +GAMMA*(SSS(I,J,IBLOCK)-ATB (I,J,1,2,IBLOCK))*seaice(i,j,IBLOCK)/ODZP(1)&
                   +GAMMA*30./365./4.*50.*(SSS(I,J,IBLOCK)-ATB (I,J,1,2,IBLOCK))/ODZP(1)*(1.0-seaice(i,j,IBLOCK)))
#else
!                     STF (I,J,IBLOCK) = GAMMA * (SSS (I,J,IBLOCK) - ATB (I,J,1,2,IBLOCK))
                      STF (I,J,IBLOCK) = GAMMA * (SSS (I,J,IBLOCK) - ATB (I,J,1,2,IBLOCK))/ODZP(1)
#endif
#endif
!                     TF (I,J,1,IBLOCK) = TF (I,J,1,IBLOCK) + STF (I,J,IBLOCK)* (1.0- AIDIF)
                     TF (I,J,1,IBLOCK) = TF (I,J,1,IBLOCK) + STF (I,J,IBLOCK)* (1.0- AIDIF)*ODZP(1)
!
!                     NET (I,J,2,IBLOCK) = STF (I,J,IBLOCK)*(1.0-AIDIF)
                     NET (I,J,2,IBLOCK) = STF (I,J,IBLOCK)*ODZP(1)
!
                  END IF
               END DO
            END DO
      END DO
          !if(mytid==100) write(16,*) 'TF,0-1',TF(21,10,1,1)
!
!LPF20160504
!===================================
!lhl20130913


#ifdef SSSNORM
! normalization for fresh water flux
       allocate(temp11(imt,jmt,max_blocks_clinic)) 
      if(.False.) allocate(temp12(imt_global,jmt_global))
      ERR_norm2=0.0
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
     DO J=3,JMT-2
     DO I=3,IMT-2
    !  DO J = JST,JMT
    !     DO I = 1,IMT
          ERR_norm2=ERR_norm2+tarea(i,j,iblock)*NET(I,J,2,iblock)*VIT(I,J,1,iblock)
     END DO
     END DO
    END DO
!$OMP END PARALLEL DO
#ifdef SPMD
      call mpi_reduce(ERR_norm2,ERR_norm1,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      !call mpi_reduce(temp11,ERR_norm1,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)

      if(.False.) then
      call gather_global(temp12, temp11, master_task, distrb_clinic) 
      if(mytid==0) then
        ERR_norm1=0.D0
!$OMP PARALLEL DO PRIVATE (J,I)
       do j=1,jmt_global
          do i=1,imt_global
            ERR_norm1=ERR_norm1+temp12(i,j)
          enddo
       enddo
     endif
    endif !noworkforSUPERHIGH

       CALL MPI_BCAST(ERR_norm1,1,MPI_PR,0,MPI_COMM_OCN,IERR)
        ERR_norm2 = - ERR_norm1/area_t

!      if (mytid.eq.0) then
!      write(*,'(a60,2D25.15)') "Gross of fresh water flux (total psu*m/s,
!      averaged psu/s) =",err1,err2
!      endif


#else
!      ERR_norm1=0.D0
!       do j=jstart+1,jmt_global-1
!          do i=2,imt_global-1
!            ERR_norm1=ERR_norm1+temp11(i,j)
!          enddo
!       enddo
!      ERR_norm2 = - ERR_norm1 / ASEA
!      write(*,'(a60,2D25.15)') "Gross of fresh water flux (total psu*m/s,
!      averaged psu/s) =",err,err2
#endif
       FW_norm2=ERR_norm2 !LPF20160823
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JMT
         DO I = 1,IMT
            TF (I,J,1,IBLOCK)=TF(I,J,1,IBLOCK)+ERR_norm2*VIT(I,J,1,IBLOCK)
         END DO
       END DO
      END DO
!$OMP END PARALLEL DO
         
         ! !if(mytid==100) then
         !   DO IBLOCK = 1, NBLOCKS_CLINIC
         !    DO J = JST,JMT
         !     DO I = 1,IMT         
         !       if(abs(tf(i,j,k,iblock))>0)write(16,*) 'TF,1-0',TF(i,j,1,iblock),i,j,iblock
         !     ENDDO
         !     ENDDO
         !   ENDDO
         !endif
          !if(mytid==100) write(16,*) 'TF,1-0',TF(21,10,1,1)

!lhl20130913
!===================================
       deallocate(temp11) 
      if(.False.) deallocate(temp12) 
#endif
!LPF20160504

  ELSE !for temperature
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
        DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = 1, JMT
               DO I = 1, IMT
                IF (KMT(I,J,IBLOCK) > 0)THEN
#ifdef COUP
                 STF (I,J,IBLOCK) = TSF(I,J,IBLOCK)
#else
#ifdef FRC_CORE
                 STF (I,J,IBLOCK) = ((SWV(I,J,IBLOCK)+NSWV(I,J,IBLOCK))*OD0CP+ & 
                                      SEAICE(I,J,IBLOCK)*GAMMA*(SST(I,J,IBLOCK)-ATB(I,J,1,1,IBLOCK))/ODZP(1))
#else
                 STF (I,J,IBLOCK) = (SWV(I,J,IBLOCK)+NSWV(I,J,IBLOCK)-DQDT(I,J,IBLOCK)* & 
                                    (SST (I,J,IBLOCK) - ATB (I,J,1,1,IBLOCK)))*OD0CP
#endif
#endif
                 TF (I,J,1,IBLOCK) = TF (I,J,1,IBLOCK) + STF (I,J,IBLOCK)* ODZP (1)* (1.0D0- AIDIF)
!
!                 NET (I,J,1,IBLOCK) = STF (I,J,IBLOCK)*ODZP(1)*(1.0-AIDIF)
                 NET (I,J,1,IBLOCK) = STF (I,J,IBLOCK)*ODZP(1)
!
                END IF !KMT
               END DO !I
            END DO !J
        END DO !IBLOCK
    END IF !N==2 
          !if(mytid==100) write(16,*) 'TF,1-1',TF(21,10,1,1)
!LPF20170729 
!for test the boudary net
#ifdef LICOM_ENABLE_TEST_TRACER
      call fortran_test_time_start ("tracer halo")
#endif
       call POP_HaloUpdate(net , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
#ifdef LICOM_ENABLE_TEST_TRACER
      call fortran_test_time_stop ("tracer halo")
#endif
!for test the boudary net
     if ( simple_assm ) then
     !if ( simple_assm==.true.) then
          if(mytid==0) write(16,*) 'into restoring,n=',N
!-----------------------------------------------------------------------
!   restoring to obs temp and salt
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
       DO K = 1,kvt
         DO j=3,jmt-2
            DO I = 3,imt-2
               TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (restore_at (I,J,K,&
                          N,iblock) - ATB (I,J,K,N,iblock))*LAMDA1(k)
                dt_restore(i,j,k,N,iblock)=VIT (I,J,K,iblock)* (restore_at (I,J,K,&
                          N,iblock) - ATB (I,J,K,N,iblock))*LAMDA1(k) 
               END DO
            END DO
         END DO
      END DO
   end if !sim_assm=true 
!LPF20170729 
!-----------------------------------------------------------------------


if ( boundary_restore == 1) then
!-----------------------------------------------------------------------
!    boundary condition
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 2,KM
         do j=3,jmt-2
               DO I = 3,imt-2
               !DO I = 3,jmt-2
                  TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (RESTORE (I,J,K,&
                             N,iblock) - ATB (I,J,K,N,iblock))*LAMDA(i,j,k,iblock)
               END DO
            END DO
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 2,KM
         do j=3,jmt-2
               DO I = 3, imt-2
                  TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (RESTORE (I,J,K,&
                             N,iblock) - ATB (I,J,K,N,iblock))*LAMDA(i,j,k,iblock)
               END DO
            END DO
         END DO
    END DO
   end if !boundary  

!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
      
          !if(mytid==100) write(16,*) 'TF,2-0',TF(21,10,1,1)
    if (adv_tracer(1:5) == 'tspas') then
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 3, JMT-2
               DO I = 3, IMT-2
                     VTL (I,J,K,iblock) = AT (I,J,K,N,iblock) + DTS * TF (I,J,K,iblock)
               END DO
            END DO
         END DO
     END DO
   else if (adv_tracer(1:8) == 'centered') then
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 3, JMT-2
               DO I = 3, IMT-2
                     VTL (I,J,K,iblock) = ATB (I,J,K,N,iblock) + C2DTTS * TF (I,J,K,iblock)
               END DO
            END DO
         END DO
     END DO
   else 
     call exit_licom(sigAbort,'The false advection option for tracer')
   end if
 
!-----------------------------------------------------------------------
!     ADD DT/DT COMPONENT DUE TO IMPLICIT VERTICAL DIFFUSION
!-----------------------------------------------------------------------

         VTL_ori=VTL !for output dt diffusion !LPF20160823
!          !if(mytid==100) write(16,*) 'TF,2-1',TF(21,10,1,1)

         CALL INVTRIT(VTL,STF,WKC,AIDIF,C2DTTS)


#ifdef LICOM_ENABLE_TEST_TRACER
      call fortran_test_time_start ("tracer halo")
#endif
     call POP_HaloUpdate(VTL , POP_haloClinic, POP_gridHorzLocCenter,&
                         POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
#ifdef LICOM_ENABLE_TEST_TRACER
      call fortran_test_time_stop ("tracer halo")
#endif
     DO K=1,KM
      if (K==1) then
       dt_diff(:,:,K,N,:)=(VTL(:,:,K,:)-VTL_ori(:,:,K,:))/C2DTTS*vit(:,:,K,:)-STF (:,:,:)*AIDIF*ODZP(1) !for output dt diffusion !LPF20160823
      else
       dt_diff(:,:,K,N,:)=(VTL(:,:,K,:)-VTL_ori(:,:,K,:))/C2DTTS*vit(:,:,K,:) 
      endif
     ENDDO
     !dt_diff(:,:,:,N,:)=(VTL(:,:,:,:)-VTL_ori(:,:,:,:))/C2DTTS*vit(:,:,:,:) !for output dt diffusion !LPF20160823
!
!---------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
 
!          !if(mytid==100) write(16,*) 'TF,3-0',TF(21,10,1,1)
               
    if (horiz_grid_opt(1:7) == 'lat_lon') then          
   
     if (mod(ist,180) == 1) then
       CALL SMTS (VTL,VIT,fil_lat2)
     else
!
      if (adv_tracer(1:5) == 'tspas') then
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
          DO IBLOCK = 1, NBLOCKS_CLINIC
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,1,IBLOCK) = VTL(I,J,1,IBLOCK)-AT (I,J,1,N,IBLOCK) - NET(I,J,N,IBLOCK)*DTS
          END DO
          END DO
          END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
          DO  IBLOCK = 1, NBLOCKS_CLINIC
          DO K=2,KM
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,K,IBLOCK) = VTL(I,J,K,IBLOCK)-AT (I,J,K,N,IBLOCK)
          END DO
          END DO
          END DO
          END DO
      else if (adv_tracer(1:8) == 'centered') then
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
          DO IBLOCK = 1, NBLOCKS_CLINIC
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,1,IBLOCK) = VTL(I,J,1,IBLOCK)- ATB (I,J,1,N,IBLOCK) - NET(I,J,N,IBLOCK)*C2DTTS
          END DO
          END DO
          END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
          DO  IBLOCK = 1, NBLOCKS_CLINIC
          DO K=2,KM
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,K,IBLOCK) = VTL(I,J,K,IBLOCK)- ATB (I,J,K,N,IBLOCK)
          END DO
          END DO
          END DO
          END DO
!
      else 
         call exit_licom(sigAbort,'The false advection option for tracer')
      end if !choose the advection scheme
!
      CALL SMTS (VTL,VIT,fil_lat1)
!
      if (adv_tracer(1:5) == 'tspas') then
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
          DO IBLOCK = 1, NBLOCKS_CLINIC
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,1,IBLOCK) = AT (I,J,1,N,IBLOCK) + VTL(I,J,1,IBLOCK) + NET(I,J,N,IBLOCK)*DTS
          END DO
          END DO
          END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
          DO  IBLOCK = 1, NBLOCKS_CLINIC
          DO K=2,KM
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,K,IBLOCK) = AT (I,J,K,N,IBLOCK) + VTL(I,J,K,IBLOCK)
          END DO
          END DO
          END DO
          END DO
      else if (adv_tracer(1:8) == 'centered') then
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
          DO IBLOCK = 1, NBLOCKS_CLINIC
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,1,IBLOCK) = ATB (I,J,1,N,IBLOCK) + VTL(I,J,1,IBLOCK) + NET(I,J,N,IBLOCK)*C2DTTS
          END DO
          END DO
          END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
          DO  IBLOCK = 1, NBLOCKS_CLINIC
          DO K=2,KM
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,K,IBLOCK) = ATB (I,J,K,N,IBLOCK) + VTL(I,J,K,IBLOCK)
          END DO
          END DO
          END DO
          END DO
      else 
         call exit_licom(sigAbort,'The false advection option for tracer')
      end if !advection scheme
   end if  !ist different
   end if !for lon-lat only 

          !if(mytid==100) write(16,*) 'TF,3-1',TF(21,10,1,1),N
!
!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1, JMT
               DO I = 1,IMT
                 tend(i,j,k,N,iblock)=(VTL(I,J,K,iblock)-AT(I,J,K,N,iblock))/C2DTTS*VIT(I,J,K,iblock)
                 !AT or ATB
               END DO
            END DO
         END DO
     END DO
          !if(mytid==100) write(16,*) 'TF,3-2,N=',TF(21,10,1,1),N

  if (adv_tracer(1:5) == 'tspas') then
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1, JMT
               DO I = 1,IMT
                  AT (I,J,K,N,IBLOCK) = VTL (I,J,K,IBLOCK)
                  ATB(I,J,K,N,IBLOCK)=AT (I,J,K,N,IBLOCK)
               END DO
            END DO
         END DO
      END DO
  else if (adv_tracer(1:8) == 'centered') then
         IF (IST >= 1)THEN
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
            DO K = 1,KM
               DO J = 1, JMT
                  DO I = 1,IMT
                     ATB (I,J,K,N,IBLOCK) = AFT2* AT (I,J,K,N,IBLOCK) + AFT1* (ATB (I,&
                                    J,K,N,IBLOCK) + VTL (I,J, K,IBLOCK))
                  END DO
               END DO
            END DO
     END DO
         END IF
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1, JMT
               DO I = 1,IMT
                  AT (I,J,K,N,IBLOCK) = VTL (I,J,K,IBLOCK)
               END DO
            END DO
         END DO
      END DO
  else 
      call exit_licom(sigAbort,'The false advection option for tracer')
  end if
  END DO
!================================================================
       
!      CALL CONVADJ !LPF20160829 !move from mainloop to here
!      DO N = 1, NTRA
!      DO K = 1,KM
!       dt_conv(:,:,k,N,:)=(AT(:,:,k,N,:)-ATB(:,:,k,N,:))/C2DTTS*vit(:,:,k,:) !for output dt diffusion !LPF20160823
!      ENDDO
!      ENDDO 
!      if (trim(adv_tracer) == 'tspas') then
!         ATB(:,:,1:km,:,:)=AT(:,:,1:km,:,:)    
!       endif
!       tend=tend+dt_conv !total tendency

!================================================================
#else
         atb(:,:,:,1,:)= 12.0D0
         atb(:,:,:,2,:)= c0
         at (:,:,:,:,:)= atb(:,:,1:km,:,:)
#endif
      IST = IST +1
!for test the boudary net
!       call POP_HaloUpdate(net , POP_haloClinic, POP_gridHorzLocCenter,&
!                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!for test the boudary net
!
#ifdef ISO
      deallocate(K1,K2,K3,adv_vetiso)
#endif
      deallocate(rict_ref,rict_replace)

      deallocate(stf,tf)
      deallocate(wkb,wkc,wkd)
      deallocate(rit)
      RETURN
      END SUBROUTINE TRACER
 
