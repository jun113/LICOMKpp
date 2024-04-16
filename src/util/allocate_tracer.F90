#include <def-undef.h>
SUBROUTINE ALLOCATE_TRACER

use work_mod
use isopyc_mod

   IMPLICIT NONE
   if (.not. allocated(stf)) then
      allocate(stf(imt, jmt, max_blocks_clinic))
   end if
   if (.not. allocated(tf)) then
      allocate(tf(imt, jmt, km, max_blocks_clinic))
   end if
   if (.not. allocated(wkb)) then
      allocate(wkb(imt, jmt, km, max_blocks_clinic))
   end if
   if (.not. allocated(wkc)) then
      allocate(wkc(imt, jmt, km, max_blocks_clinic))
   end if
   if (.not. allocated(wkd)) then
      allocate(wkd(imt, jmt, km, max_blocks_clinic))
   end if

#ifdef ISO
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
   if (.not. allocated(adv_vetiso)) then
      allocate(adv_vetiso(imt, km, jmt, max_blocks_clinic))
   end if

   ! isoflux
   if (.not. allocated(work_1)) then
      allocate(work_1(imt, jmt, km, max_blocks_clinic))
   endif
   if (.not. allocated(work_2)) then
      allocate(work_2(imt, jmt, km, max_blocks_clinic))
   endif
   if (.not. allocated(work_3)) then
      allocate(work_3(imt, jmt, 0:km, max_blocks_clinic))
   endif
   if (.not. allocated(temp)) then
      allocate(temp(imt, jmt, km, max_blocks_clinic))
   endif
#endif 
! ifdef ISO

END SUBROUTINE ALLOCATE_TRACER
