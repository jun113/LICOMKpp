SUBROUTINE ALLOCATE_READYT

use dyn_mod
use pmix_mod

   IMPLICIT NONE

   if (.not. allocated(dlu)) then
      allocate(dlu(imt, jmt, km, max_blocks_clinic))
   end if
   if (.not. allocated(dlv)) then
      allocate(dlv(imt, jmt, km, max_blocks_clinic))
   end if
   if (.not. allocated(gg)) then
      allocate(gg(imt, jmt, km, max_blocks_clinic))
   end if
   if (.not. allocated(ric)) then
      allocate(ric(imt, jmt, kmm1, max_blocks_clinic))
   end if
   if (.not. allocated(rit)) then
      allocate(rit(imt, jmt, kmm1, max_blocks_clinic))
   end if
   if (.not. allocated(rict)) then
      allocate(rict(imt, jmt, kmm1, max_blocks_clinic))
   end if
   if (.not. allocated(rict_replace)) then
      allocate(rict_replace(imt, jmt, kmm1, max_blocks_clinic))
   end if

END SUBROUTINE
