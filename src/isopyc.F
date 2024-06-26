!  CVS: $Id: isopyc.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE ISOPYC
!     =================
 
!     Compute the isopycnal mixing tensor components and the
!     isopycnal advection velocities which param_modeterize the effect
!     of eddies on the isopycnals.
 
!     Mixing tensor "K" is ...
!          | 1.0            K1(,,,2)          K1(,,,3) |
!          |                                           |
!     K =  | K2(,,,1)        1.0              K2(,,,3) |
!          |                                           |
!          | K3(,,,1)       K3(,,,2)          K3(,,,3) |
 
!     where K1(,,,2) and K2(,,,1) are set to 0.0 (neglected)
 
!     output:
!       rhoi = density at tau-1 referenced to pressure levels
!       K1   = tensor components (1,2), (1,3) centered on east face
!              of "T" cells
!       K2   = tensor components (2,1), (2,3) centered on north face
!              of "T" cells
!       K3   = tensor components (3,1), (3,2), (3,3) centered on
!              bottom face of "T" cells
!       adv_vetiso = isopycnal advective vel on east face of "T" cell
!       adv_vntiso = isopycnal advective vel on north face of "T" cell
!               (Note: this includes the cosine factor as in "adv_vnt")
!       adv_vbtiso = isopycnal advective vel on bottom face of "T" cell
use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use isopyc_mod
use domain
      IMPLICIT NONE
 
      REAL(r8):: tq,sq,tref0,sref0
      INTEGER :: kq, iblock
      REAL(r8):: DENS
      EXTERNAL DENS
 
!-----------------------------------------------------------------------
!     compute normalized densities for each isopycnal reference pressure
!     level using a 3rd order polynomial fit to the equation of state.
!     for each isopycnal reference pressure level, the same reference
!     potential temperature, reference salinity and expansion coeff
!     values are used at all of the vertical levels.
 
!     Note: this density is used for the mixing tensor in both the
!     Redi/Cox and Gent/McWilliams options
!-----------------------------------------------------------------------
!     wjl 20211031 
   if (.not. allocated(e)) then
      allocate(e(imt, kmp1, jmt, 3, max_blocks_clinic))
   end if
   if (.not. allocated(rhoi)) then
      allocate(rhoi(imt, 0:km, jmt, nrpl, max_blocks_clinic))
   end if
   if (.not. allocated(k1)) then
      allocate(k1(imt, 0:km, jmt, 3:3, max_blocks_clinic))
   end if
   if (.not. allocated(k2)) then
      allocate(k2(imt, 0:km, jmt, 3:3, max_blocks_clinic))
   end if
   if (.not. allocated(k3)) then
      allocate(k3(imt, 0:km, jmt, 1:3, max_blocks_clinic))
   end if
!      allocate(e(imt,kmp1,jmt,3,max_blocks_clinic),rhoi(imt,0:km,jmt,nrpl,max_blocks_clinic))
!      allocate(K1(imt,0:km,jmt,3:3,max_blocks_clinic),K2(imt,0:km,jmt,3:3,max_blocks_clinic), & 
!               K3(imt,0:km,jmt,1:3,max_blocks_clinic))


    rhoi(:,0,:,:,:) = 0.0D0


!$OMP PARALLEL DO PRIVATE (iblock,m,tref0,sref0,tq,sq,kq)
  do iblock = 1, nblocks_clinic
      DO m = 1,nrpl
         tref0 = to (krplin (m))
         sref0 = so (krplin (m))
       DO k = 1,km
          DO j = 1,jmt
               DO i = 1,imt
                  tq = atb (i,j,k,1,iblock) - tref0
                  sq = atb (i,j,k,2,iblock) - sref0
                  kq = krplin (m)
                  ! wjl 20211103
                  rhoi (i,k,j,m,iblock) = dens (tq, sq, kq)
!                 rhoi (i,k,j,m,iblock) = dens (tq, sq, kq,iblock)
               END DO
            END DO
         END DO
      END DO
  end do

 
!-----------------------------------------------------------------------
!     evaluate K2(,,3) centered on the northern face of "T" cells
!-----------------------------------------------------------------------
 
      CALL k2_3
 
!-----------------------------------------------------------------------
!     evaluate K1(,,3) centered on eastern face of "T" cells
!-----------------------------------------------------------------------
 
      CALL k1_3
 
!-----------------------------------------------------------------------
!     evaluate K3(,,1..3) centered on bottom face of "T" cells
!-----------------------------------------------------------------------
 
      CALL k3_123
 
!-----------------------------------------------------------------------
!     compute isopycnal advective velocities for tracers
!-----------------------------------------------------------------------
 
      deallocate(e,rhoi)
! wjl 20211101
      if (.not. allocated(adv_vetiso)) then
         allocate(adv_vetiso(imt,km,jmt,max_blocks_clinic) )
      end if

!      allocate(adv_vetiso(imt,km,jmt,max_blocks_clinic) )

      CALL isoadv
 
      RETURN
      END SUBROUTINE ISOPYC
 
 
#else
      SUBROUTINE ISOPYC ()
      RETURN
 
      END SUBROUTINE ISOPYC
 
#endif 
