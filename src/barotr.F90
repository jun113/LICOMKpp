      SUBROUTINE BAROTR
!     =================
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use msg_mod
use domain
use grid
use blocks
use hmix_del2
use hmix_del4
use operators
use smuvh
use POP_GridHorzMod
use POP_HaloMod
use global_reductions
use gather_scatter
use distribution
use constant_mod
      IMPLICIT NONE

      INTEGER :: IEB,NC,IEB_LOOP
      real(r8)    :: gstar ,am_viv,fil_lat1,fil_lat2, ek0, maxz0, minz0
!     real(r8)  ::  ttt(imt_global,jmt_global)
      integer :: iblock,ii1,jj1,ii2,jj2, irec
      real(r8):: hduk(imt,jmt) , hdvk(imt,jmt), gradx(imt,jmt),grady(imt,jmt), div_out(imt,jmt)
      real(r8):: hdtk(imt,jmt)
      REAL (r8)   :: DT2K(imt,jmt)
      type(block):: this_block

!      Define the threthold latitute for zonal smoother
      allocate (buffer(imt_global,jmt_global))
       fil_lat1=65.0D0
       fil_lat2=65.0D0
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
#ifdef LICOM_ENABLE_TEST_TIME
#undef  LICOM_ENABLE_TEST_BAROTR
#define LICOM_ENABLE_TEST_BAROTR
#endif
      wka=0
      work=0
      baro_loop : DO NC = 1,NBB

!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
        call tgrid_to_ugrid(work(:,:,iblock),h0(:,:,iblock),iblock)
   end do
!---------------------------------------------------------------------
!     COMPUTING DH0
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (iblock,J)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 2, JMT
            DO I = 1,IMT-1
               WKA (I,J,1,IBLOCK)= UB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
               WKA (I,J,2,IBLOCK)= VB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
            END DO 
         END DO
     END DO
     

!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,J,div_out) 
    DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call div(1,DIV_OUT,wka(:,:,1,iblock),wka(:,:,2,iblock),this_block)
         DO J = 3, jmt-2
            DO I = 3,imt-2
               WORK (I,J,IBLOCK)=VIT(I,J,1,IBLOCK)*(-1)*div_out(i,j)*P25
             END DO 
          ENDDO
    END DO
!
#ifdef LICOM_ENABLE_TEST_BAROTR
      !call fortran_test_time_start ("barotr halo work")
      call fortran_test_time_start ("barotr halo")
#endif
         call POP_HaloUpdate(work , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

#ifdef LICOM_ENABLE_TEST_BAROTR
      !call fortran_test_time_stop ("barotr halo work")
      call fortran_test_time_stop ("barotr halo")
#endif
    
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1, jmt
            DO I = 1,imt
               H0 (I,J,IBLOCK)= H0P (I,J,IBLOCK) + WORK (I,J,IBLOCK) * DTB
            END DO
         END DO
     END DO
!
!---------------------------------------------------------------------
!     COMPUTE THE "ARTIFICIAL" HORIZONTAL VISCOSITY
!---------------------------------------------------------------------

#if ( defined SMAG1)
#else
#if (defined BIHAR)
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block, hduk,hdvk)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call hdiffu_del4(1, HDUK, HDVK, ubp(:,:,iblock), vbp(:,:,iblock), this_block)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,5,IBLOCK)= hduk(i,j)
               WKA (I,J,6,IBLOCK)= hdvk(i,j)
            END DO
         END DO
     END DO
!
!!!!!!
#else
!$OMP PARALLEL DO PRIVATE (IBLOCK, this_block, hduk, hdvk)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call hdiffu_del2(1, HDUK, HDVK, ubp(:,:,iblock), vbp(:,:,iblock), this_block)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,5,IBLOCK)= hduk(i,j)
               WKA (I,J,6,IBLOCK)= hdvk(i,j)
            END DO
         END DO
    END DO

#endif
#endif
!      if (mytid == 97) then
!         buffer_real4(:,:)= hduk(:,:)
!         irec=(NC-1)*8+1
!         write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!      end if
!
!     + (g'-1)g*dH/dr
!---------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (IBLOCK, this_block, gradx, grady, gstar)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call grad(1, GRADX, GRADY, H0(:,:,iblock), this_block)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               gstar=(WGP (I,J,IBLOCK) -1.0)*G 
               WKA (I,J,1,IBLOCK) = WKA (I,J,5,IBLOCK) + gstar*GRADX(I,J)
               WKA (I,J,2,IBLOCK) = WKA (I,J,6,IBLOCK) + gstar*GRADY(I,J)
            END DO
         END DO
     END DO

!Yu
!      if (mytid == 97) then
!         buffer_real4(:,:)= wka(:,:,1,1)
!         irec=(NC-1)*8+2
!         write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!      end if
!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
        call tgrid_to_ugrid(work(:,:,iblock),h0(:,:,iblock),iblock)
   end do
!
!---------------------------------------------------------------------
!     COMPUTING DU & DV
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,J)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,1,IBLOCK)= VIV (I,J,1,IBLOCK)* ( WKA (I,J,1,IBLOCK) + DLUB (I,J,IBLOCK)     &
                              - FCOR(I,J,Iblock)* VBP (I,J,IBLOCK) + &
               PAX (I,J,IBLOCK) + PXB (I,J,IBLOCK) - WORK (I,J,IBLOCK)* WHX (I,J,IBLOCK) )
               WKA (I,J,2,IBLOCK)= VIV (I,J,1,IBLOCK)* ( WKA (I,J,2,IBLOCK) + DLVB (I,J,IBLOCK)     &
                              + FCOR(I,J,Iblock)* UBP (I,J,IBLOCK) + &
               PAY (I,J,IBLOCK) + PYB (I,J,IBLOCK) - WORK (I,J,IBLOCK)* WHY (I,J,IBLOCK) )
            END DO
         END DO
     END DO

!---------------------------------------------------------------------
!     CORIOLIS ADJUSTMENT
!---------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (IBLOCK,J)
     DO IBLOCK = 1, NBLOCKS_CLINIC
     DO J = 3, jmt-2
     DO I = 3, imt-2
        WKA (I,J,3,IBLOCK)= EBEA(I,J,IBLOCK)* WKA (I,J,1,IBLOCK) - EBEB(I,J,IBLOCK)* WKA (I,J,2,IBLOCK)
        WKA (I,J,4,IBLOCK)= EBEA(I,J,IBLOCK)* WKA (I,J,2,IBLOCK) + EBEB(I,J,IBLOCK)* WKA (I,J,1,IBLOCK)
     END DO
     END DO
     END DO
!
#ifdef LICOM_ENABLE_TEST_BAROTR
      !call fortran_test_time_start ("barotr halo wka")
      call fortran_test_time_start ("barotr halo")
#endif
    call POP_HaloUpdate(wka(:,:,3,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                        POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
    call POP_HaloUpdate(wka(:,:,4,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                        POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
#ifdef LICOM_ENABLE_TEST_BAROTR
      !call fortran_test_time_stop ("barotr halo wka")
      call fortran_test_time_stop ("barotr halo")
#endif

!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1, jmt
            DO I = 1,imt
               UB (I,J,IBLOCK)= UBP (I,J,IBLOCK) + WKA (I,J,3,IBLOCK)* DTB
               VB (I,J,IBLOCK)= VBP (I,J,IBLOCK) + WKA (I,J,4,IBLOCK)* DTB
            END DO
         END DO
     END DO
!---------------------------------------------------------------------
!     COMPUTING DH0
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 2, JMT
            DO I = 1,IMT-1
               WKA (I,J,1,IBLOCK)= UB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
               WKA (I,J,2,IBLOCK)= VB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
            END DO
         END DO
     END DO

!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,hdtk,div_out) 
    DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call div(1,DIV_OUT,wka(:,:,1,iblock),wka(:,:,2,iblock),this_block)
!
!2020 YYQ
         !if ( mod(nc,2) == 0 ) then
         if ( mod(nc,4) == 0 ) then
!2020 YYQ
#ifdef  BIHAR
!2020 YYQ
             !call hdifft_del4(1,hdtk,h0p(:,:,iblock),this_block)
             call hdifft_del4(1,dt2k,hdtk,h0p(:,:,iblock),this_block)
!2020 YYQ
#else
             call hdifft_del2(1,hdtk,h0p(:,:,iblock),this_block)
#endif
         else
             hdtk = c0
         end if
!
         DO J = 2, jmt-2
            DO I = 3,imt-2
               !WORK (I,J,IBLOCK)=VIT(I,J,1,IBLOCK)*(hdtk(i,j) - div_out(i,j))
               WORK (I,J,IBLOCK)=VIT(I,J,1,IBLOCK)*(hdtk(i,j)*1.0D0 - div_out(i,j))
             END DO
          ENDDO
    END DO
!
!---------------------------------------------------------------------
!     PREDICTING VB , UB & H0
!---------------------------------------------------------------------
#ifdef LICOM_ENABLE_TEST_BAROTR
      !call fortran_test_time_start ("barotr halo work")
      call fortran_test_time_start ("barotr halo")
#endif
         call POP_HaloUpdate(work , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
#ifdef LICOM_ENABLE_TEST_BAROTR
      !call fortran_test_time_stop ("barotr halo work")
      call fortran_test_time_stop ("barotr halo")
#endif
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1, jmt
            DO I = 1,imt
               H0 (I,J,IBLOCK)= H0P (I,J,IBLOCK) + WORK (I,J,IBLOCK) * DTB
            END DO
         END DO
     END DO
 
            ISB = ISB +1
            ubp=ub
            vbp=vb
            h0p=h0
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = JST,JET ! Dec. 4, 2002, Yongqiang YU
               DO I = 1,IMT
                  H0F (I,J,IBLOCK) = H0F (I,J,IBLOCK) + H0 (I,J,IBLOCK)
                  H0BF (I,J,IBLOCK) = H0BF (I,J,IBLOCK) + H0 (I,J,IBLOCK)
               END DO
            END DO
    END DO
!    if (maxval(abs(ub)) > 5.0D0 ) then
!         write(6,*) "UBMAX", isc, isb, mytid,maxval(abs(ub)),maxloc(abs(ub))
!         write(6,'(1x,A5,1x,3I8,1x,F24.20,1x,4I4)') "UBMAX", isc, isb, mytid,maxval(abs(ub)),maxloc(abs(ub))
!ZWP20170209
!    end if
!    if (maxval(abs(vb)) > 5.0D0 ) then
!         write(6,*) "VBMAX", isc, isb, mytid,maxval(abs(vb)),maxloc(abs(vb))
!         write(6,'(1x,A5,1x,3I8,1x,F24.20,1x,4I4)') "UBMAX", isc, isb, mytid,maxval(abs(ub)),maxloc(abs(ub))
!ZWP20170209
!    end if

      END DO baro_loop
    deallocate (buffer)

      deallocate(dlub,dlvb)
      RETURN
      END SUBROUTINE BAROTR


