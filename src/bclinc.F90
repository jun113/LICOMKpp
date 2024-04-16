!  CVS: $Id: bclinc.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE BCLINC
!     =================
!     INTEGRATION OF MOMENTUM EQUATION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use forc_mod, only: psa,su,sv
use domain
use grid
use blocks
use smuvh
use operators

      IMPLICIT NONE
      REAL(r8)    :: AA,GGU,WK1,WK2,fil_lat1,fil_lat2,eee
      integer     :: iblock,kt,irec,kmb
      REAL(r8)    :: gradx(imt,jmt), grady(imt,jmt)
      REAL(r8)    :: ucheck(imt,jmt,km), vcheck(imt,jmt,km)!ZWP20170212
      integer     :: zwploc(3)
      type (block) :: this_block

#ifdef LICOM_ENABLE_TEST_TIME
#undef  LICOM_ENABLE_TEST_BCLINC
#define LICOM_ENABLE_TEST_BCLINC
#endif

!lhl0711
#ifdef CANUTO
      REAL(r8)    :: AIDIF

!M
 
 
!$OMP PARALLEL DO PRIVATE (iblock,J,I) 
      do iblock = 1, nblocks_clinic
         DO J = 2, jmt-1
            DO I = 2,imt-1
              kmb = kmu(i,j,iblock)
              if( kmb>=1 ) then !LPF20170621
               SBCX(I,J,IBLOCK) = SU (I,J,IBLOCK)* OD0
               SBCY(I,J,IBLOCK) = SV (I,J,IBLOCK)* OD0
               BBCX(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))&
                          *(UP(I,J,kmb,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP (I,J,kmb,IBLOCK)*SAG)
               BBCY(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))&
                          *(-SNLAT(I,J,IBLOCK)*UP(I,J,kmb,IBLOCK)*SAG+VP(I,J,kmb,IBLOCK)*CAG)
              end if
            ENDDO
         ENDDO
     ENDDO

      AIDIF=0.5  
!      if (ISC/=0)  AIDIF=0.5
 
#endif
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
       fil_lat1=63.0D0
       fil_lat2=63.0D0
!lhl0711

!---------------------------------------------------------------------
!     ADVECTION + DIFFUSION + CORIOLIS
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 2,jmt-1
            DO I = 2,imt-1
               DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) - FCOR(I,J,IBLOCK) * VP (I,J,K,IBLOCK)
               DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + FCOR(I,J,IBLOCK) * UP (I,J,K,IBLOCK)
!              buffer_real4(i,j)= - FCOR(I,J,IBLOCK) * VP (I,J,1,IBLOCK)
            END DO
         END DO
      END DO
  END DO
!    if (mytid ==68 ) then
!       irec = 16*isc +2
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!    end if
!
!---------------------------------------------------------------------
!     PRESSURE GRADIENT FORCES
!---------------------------------------------------------------------
!  if (isc < 10 .and. mytid == 5) then
!    write(100+mytid,*) "OK--2, ISC=", ISC, dlu(14,16,1,1), dlv(14,16,1,1)
!  end if
 
!!@@@@ DP'/DX
 
      AA = 0.0
      IF (ISC /= 0) AA = 0.5
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1,jmt
         DO I = 1,IMT
            H0BF (I,J,IBLOCK)= H0BF (I,J,IBLOCK)* ONBB
         END DO
      END DO
  END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1, jmt
         DO I = 1,IMT
            WORK (I,J,IBLOCK)= AA * H0BF (I,J,IBLOCK) + (1.0- AA)* H0BL (I,J,IBLOCK)
         END DO
      END DO
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,WKK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1,jmt
         DO I = 1,IMT
            WKK (1)= (PSA (I,J,IBLOCK)* OD0+ WORK (I,J,IBLOCK)* G)* VIT (I,J,1,IBLOCK)
            DO K = 1,KM
               WKK (K +1)= WKK (K) - GG (I,J,K,IBLOCK)* DZP (K)* VIT (I,J,K,IBLOCK)
            END DO
 
            DO K = 1,KM
               WKA (I,J,K,IBLOCK)= P5* (WKK (K) + WKK (K +1))
            END DO
         END DO
      END DO
  END DO
 
 
!$OMP PARALLEL DO PRIVATE (iblock, this_block,gradx,grady,K)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call grad(k, GRADX, GRADY, wka(:,:,k,iblock), this_block)
         DO J = 2,jmt-1
            DO I = 2,imt-1
               DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) - grady(i,j)
               DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) - gradx(i,j)
!              if (k ==1 ) buffer_real4(i,j) =  -gradx(i,j)
            END DO
         END DO
      END DO
  END DO
!
!!@@@@ G'DH/DX
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= (1.0+ OHBT (I,J,IBLOCK)* ZKT (K))* WORK (I,J,IBLOCK)      &
                            * VIT (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block, K,GGU, gradx,grady)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call grad(k, GRADX, GRADY, wka(:,:,k,iblock), this_block)
         DO J = 2, JMT-1
            DO I = 2,IMT-1
               GGU = P25* (GG (I,J,K,IBLOCK) + GG (I -1,J,K,IBLOCK) + GG (I,J +1,K,IBLOCK) &
                     + GG (I -1,J +1,K,IBLOCK))
               DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + GGU * grady(i,j)
               DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) + GGU * gradx(i,j)
!              if (k==1) buffer_real4(i,j)= GGU*gradx(i,j)
            END DO
         END DO
      END DO
   END DO

!---------------------------------------------------------------------
!     CORIOLIS ADJUSTMENT
!---------------------------------------------------------------------
 
      IF (ISC == 0) THEN
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,WK1,WK2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 2, JMT-1
               DO I = 2,IMT-1
                  WK1 = EPEA(I,J,IBLOCK)* DLV(I,J,K,IBLOCK) + EPEB(I,J,IBLOCK)* DLU (I,J,K,IBLOCK)
                  WK2 = EPEA(I,J,IBLOCK)* DLU(I,J,K,IBLOCK) - EPEB(I,J,IBLOCK)* DLV (I,J,K,IBLOCK)
                  DLV (I,J,K,IBLOCK)= WK1* VIV (I,J,K,IBLOCK)
                  DLU (I,J,K,IBLOCK)= WK2* VIV (I,J,K,IBLOCK)

               END DO
            END DO
         END DO
  END DO
      ELSE
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,WK1,WK2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 2, JMT-1
               DO I = 2,IMT-1
                  WK1 = EPLA(I,J,IBLOCK)* DLV(I,J,K,IBLOCK) + EPLB(I,J,IBLOCK)*DLU(I,J,K,IBLOCK)
                  WK2 = EPLA(I,J,IBLOCK)* DLU(I,J,K,IBLOCK) - EPLB(I,J,IBLOCK)*DLV(I,J,K,IBLOCK)
                  DLV (I,J,K,IBLOCK)= WK1* VIV(I,J,K,IBLOCK)
                  DLU (I,J,K,IBLOCK)= WK2* VIV(I,J,K,IBLOCK)

               END DO
            END DO
         END DO
  END DO
      END IF

    if (trim(horiz_grid_opt) == 'lat_lon') then
         call POP_HaloUpdate(DLU, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
         call POP_HaloUpdate(DLV, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
      CALL SMUV_3D (DLV,VIV,fil_lat1)
      CALL SMUV_3D (DLV,VIV,fil_lat1)
    end if
!YU 

!---------------------------------------------------------------------
!     PREDICTING VC & UC
!---------------------------------------------------------------------
 
      IF (ISC < 1) THEN
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 3,JMT-2
            DO I = 3,IMT-2
               V (I,J,K,IBLOCK)= VP (I,J,K,IBLOCK) + DLV (I,J,K,IBLOCK)* DTC
               U (I,J,K,IBLOCK)= UP (I,J,K,IBLOCK) + DLU (I,J,K,IBLOCK)* DTC
            END DO
         END DO
      END DO
  END DO
!
!    if (mytid ==68 ) then
!       buffer_real4(:,:)= u(:,:,1,1)
!       irec = 16*isc +12
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!    end if
!    write(mytid+120,*) ((u(i,j,30,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((sbcx(i,j,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((bbcx(i,j,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((akmu(i,j,29,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((akmu(i,j,30,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((kmu(i,j,2),i=3,imt-2),j=3,jmt-2)
     CALL INVTRIU(U,SBCX,BBCX,AKMU,AIDIF,DTC)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((u(i,j,30,2),i=3,imt-2),j=3,jmt-2)
!    close(120+mytid)
     CALL INVTRIU(V,SBCY,BBCY,AKMU,AIDIF,DTC)
!    stop
!
!    if (mytid ==68 ) then
!       buffer_real4(:,:)= u(:,:,1,1)
!       irec = 16*isc +13
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!    end if
#ifdef LICOM_ENABLE_TEST_BCLINC
      call fortran_test_time_start ("bclinc halo")
#endif
         call POP_HaloUpdate(U, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
         call POP_HaloUpdate(V, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
#ifdef LICOM_ENABLE_TEST_BCLINC
      call fortran_test_time_stop ("bclinc halo")
#endif
!
!  if (isc < 10 .and. mytid == 5) then
!    write(100+mytid,*) "OK--6, ISC=", ISC, u(14,16,1,1), v(14,16,1,1)
!  end if
!
!---------------------------------------------------------------------
!@@@  INTERACTION BETWEEN BAROTROPIC AND BAROCLINIC MODES
!---------------------------------------------------------------------
      CALL VINTEG (U,WORK)
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               U (I,J,K,IBLOCK)= (U (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + UB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
!    if (mytid ==68 ) then
!       buffer_real4(:,:)= work(:,:,1)
!       irec = 16*isc +14
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!       buffer_real4(:,:)= ub(:,:,1)
!       irec = 16*isc +15
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!       buffer_real4(:,:)= gg(:,:,1,1)
!       irec = 16*isc +16
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!    end if
!
      CALL VINTEG (V,WORK)
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               V (I,J,K,IBLOCK)= (V (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + VB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO

!  if (isc < 10 .and. mytid == 5) then
!    write(100+mytid,*) "OK--6, ISC=", ISC, u(14,16,1,1), v(14,16,1,1)
!  end if
 
      ISC = ISC +1
    
      ELSE 

 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 3,JMT-2
            DO I = 3,IMT-2
               WKA (I,J,K,IBLOCK)= VP (I,J,K,IBLOCK) + DLV (I,J,K,IBLOCK)* DTC2
            END DO
         END DO
      END DO
   END DO
!
     CALL INVTRIU(WKA,SBCY,BBCY,AKMU,AIDIF,DTC2)
!
#ifdef LICOM_ENABLE_TEST_BCLINC
      call fortran_test_time_start ("bclinc halo")
#endif
         call POP_HaloUpdate(WKA, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
#ifdef LICOM_ENABLE_TEST_BCLINC
      call fortran_test_time_stop ("bclinc halo")
#endif

  
 
 
      CALL VINTEG (WKA,WORK)
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= (WKA (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + VB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
 
!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               VP (I,J,K,IBLOCK) = AFC2* V (I,J,K,IBLOCK) + AFC1* (VP (I,J,K,IBLOCK) + WKA (I,J,K,IBLOCK))
               V (I,J,K,IBLOCK) = WKA (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 3,jmt-2
            DO I = 3,imt-2
               WKA (I,J,K,IBLOCK)= UP (I,J,K,IBLOCK) + DLU (I,J,K,IBLOCK)* DTC2
            END DO
         END DO
      END DO
  END DO
!
     CALL INVTRIU(WKA,SBCX,BBCX,AKMU,AIDIF,DTC2)
!
#ifdef LICOM_ENABLE_TEST_BCLINC
      call fortran_test_time_start ("bclinc halo")
#endif
         call POP_HaloUpdate(WKA, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
#ifdef LICOM_ENABLE_TEST_BCLINC
      call fortran_test_time_stop ("bclinc halo")
#endif
      CALL VINTEG (WKA,WORK)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= (WKA (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + UB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               UP (I,J,K,IBLOCK) = AFC2* U (I,J,K,IBLOCK) + AFC1* (UP (I,J,K,IBLOCK) + WKA (I,J,K,IBLOCK))
               U (I,J,K,IBLOCK) = WKA (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
!YU  Oct. 24, 1605
!lhl0711     IF (MOD(ISC,60)==0) THEN
    if (trim(horiz_grid_opt) == 'lat_lon') then
     IF (MOD(ISC,160)==1) THEN
        CALL SMUV_3D (U,VIV,fil_lat2)
        CALL SMUV_3D (V,VIV,fil_lat2)
        CALL SMUV_3D (UP,VIV,fil_lat2)
        CALL SMUV_3D (VP,VIV,fil_lat2)
     END IF
    end if
!
      ISC = ISC +1
      END IF
!     deallocate(buffer_real4)
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) 
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UTF (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK) + U (I,J,K,IBLOCK)
               VTF (I,J,K,IBLOCK)= VTF (I,J,K,IBLOCK) + V (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO

!ZWP20170212 
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
!      write(253,*)'isc=',isc,'mytid=',mytid,'i_global=',this_block%i_glob,'j_global=',this_block%j_glob
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
      ucheck(i,j,k) = u(i,j,k,iblock) 
      vcheck(i,j,k) = v(i,j,k,iblock) 
            END DO
         END DO
      END DO
!ZWP20170212 
  if (maxval(abs(ucheck)) > 5.0D0 ) then
         zwploc = maxloc(abs(ucheck))
!        write(6,*) "UMAX", isc, mytid,maxval(abs(u)),maxloc(abs(u))
        write(6,'(1x,A5,1x,2I8,1x,F24.20,1x,3I3,1x,2I6)') "UMAX",isc, mytid,&
             maxval(abs(ucheck)),maxloc(abs(ucheck)),&
             this_block%i_glob(zwploc(1)),this_block%j_glob(zwploc(2))
!ZWP20170209
  end if
  if (maxval(abs(vcheck)) > 5.0D0 ) then
         zwploc = maxloc(abs(vcheck))
!        write(6,*) "VMAX", isc, mytid,maxval(abs(v)),maxloc(abs(v))
        write(6,'(1x,A5,1x,2I8,1x,F24.20,1x,3I3,1x,2I6)') "VMAX",isc, mytid,&
             maxval(abs(vcheck)),maxloc(abs(vcheck)),&
             this_block%i_glob(zwploc(1)),this_block%j_glob(zwploc(2))
  end if

  END DO
!ZWP2017dd0209

  call mpi_barrier(mpi_comm_ocn,ierr)

      RETURN
      END SUBROUTINE BCLINC
