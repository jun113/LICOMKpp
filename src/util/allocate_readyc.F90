SUBROUTINE ALLOCATE_READYC
use dyn_mod
use pmix_mod
use work_mod

IMPLICIT NONE
   if (.not. allocated(riu)) then
      allocate(riu(imt, jmt, 0:km, max_blocks_clinic), stat = ierr) 
   end if
   if (.not. allocated(rict_ref)) then
      allocate(rict_ref(imt, jmt, max_blocks_clinic), stat = ierr) 
   end if 
   if (.not. allocated(dlub)) then
      allocate(dlub(imt, jmt, max_blocks_clinic), stat = ierr) 
   end if 
   if (.not. allocated(dlvb)) then
      allocate(dlvb(imt, jmt, max_blocks_clinic), stat = ierr) 
   end if 

   if (.not. allocated(uk)) then
      allocate(uk(imt, jmt, km, max_blocks_clinic), stat = ierr) 
   end if 
   if (.not. allocated(vk)) then
      allocate(vk(imt, jmt, km, max_blocks_clinic), stat = ierr) 
   end if 
END SUBROUTINE ALLOCATE_READYC
