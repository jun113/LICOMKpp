subroutine scatter_global_jra_r4(val, buffer, tid)
use precision_mod,  only: i4, r4
use domain,         only: distrb_clinic
use gather_scatter, only: scatter_global
use constant_mod,   only: &
   field_loc_center, field_type_scalar
use param_mod,      only: &
   imt, jmt, max_blocks_clinic, imt_global, jmt_global, master_task
implicit none
  real(r4), dimension(imt, jmt, max_blocks_clinic),intent(inout) :: val
	real(r4), dimension(imt_global, jmt_global), intent(in) :: buffer
	integer(i4) :: tid

   call scatter_global(val, buffer, tid, distrb_clinic, &
	    field_loc_Center, field_type_scalar)

end subroutine scatter_global_jra_r4

subroutine scatter_global_jra_r8(val, buffer, tid)
use precision_mod,  only: i4, r8
use domain,         only: distrb_clinic
use gather_scatter, only: scatter_global
use constant_mod,   only: &
   field_loc_center, field_type_scalar
use param_mod,      only: &
   imt, jmt, max_blocks_clinic, imt_global, jmt_global, master_task
implicit none
  real(r8), dimension(imt, jmt, max_blocks_clinic),intent(inout) :: val
	real(r8), dimension(imt_global, jmt_global), intent(in) :: buffer
	integer(i4) :: tid

   call scatter_global(val, buffer, tid, distrb_clinic, &
	    field_loc_Center, field_type_scalar)

end subroutine scatter_global_jra_r8