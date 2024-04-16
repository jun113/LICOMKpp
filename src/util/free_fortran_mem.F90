
subroutine free_forc_mod
use forc_mod
implicit none
   if (allocated(chloro3)) then
	    deallocate(chloro3)
	 endif 
   if (allocated(su3_io)) then
	    deallocate(su3_io)
	 endif 
   if (allocated(psa3_io)) then
	    deallocate(psa3_io)
	 endif 
   if (allocated(tsa3_io)) then
	    deallocate(tsa3_io)
	 endif 
   if (allocated(qar3_io)) then
	    deallocate(qar3_io)
	 endif 
   if (allocated(uva3_io)) then
	    deallocate(uva3_io)
	 endif 
   if (allocated(swv3_io)) then
	    deallocate(swv3_io)
	 endif 
   if (allocated(cld3_io)) then
	    deallocate(cld3_io)
	 endif 
   if (allocated(sss3_io)) then
	    deallocate(sss3_io)
	 endif 
   if (allocated(sst3_io)) then
	    deallocate(sst3_io)
	 endif 
   if (allocated(nswv3_io)) then
	    deallocate(nswv3_io)
	 endif 
   if (allocated(dqdt3_io)) then
	    deallocate(dqdt3_io)
	 endif 
   if (allocated(chloro3_io)) then
	    deallocate(chloro3_io)
	 endif 
   if (allocated(wspdu3_io)) then
	    deallocate(wspdu3_io)
	 endif 
   if (allocated(wspdv3_io)) then
	    deallocate(wspdv3_io)
	 endif 
   if (allocated(lwv3_io)) then
	    deallocate(lwv3_io)
	 endif 
   if (allocated(seaice3_io)) then
	    deallocate(seaice3_io)
	 endif 
   if (allocated(rain3_io)) then
	    deallocate(rain3_io)
	 endif 
   if (allocated(snow3_io)) then
	    deallocate(snow3_io)
	 endif 
   if (allocated(runoff3_io)) then
	    deallocate(runoff3_io)
	 endif 
   if (allocated(buffer3d)) then
	    deallocate(buffer3d)
	 endif 
   if (allocated(w3d)) then
	    deallocate(w3d)
	 endif 
end subroutine free_forc_mod

subroutine free_grid
use grid
implicit none
   if (allocated(ht)) then
      deallocate(ht) 
   end if
   if (allocated(hu)) then
      deallocate(hu) 
   end if
end subroutine free_grid

subroutine free_work_mod
use work_mod
implicit none
   if (allocated(uk)) then
      deallocate(uk) 
   end if
   if (allocated(vk)) then
      deallocate(vk) 
   end if
end subroutine free_work_mod