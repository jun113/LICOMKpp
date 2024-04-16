!  CVS: $Id: dyn_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module dyn_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------
      real(r8),dimension(imt,jmt,max_blocks_clinic)::ub,vb,ubp,vbp,h0p
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::up,vp
      real(r8),dimension(imt,jmt,kmp1,max_blocks_clinic)::ws
      real(r8),dimension(imt,jmt,max_blocks_clinic)::h0l,h0f,h0bl,h0bf
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::utl,utf,vtl,vtf
!
      REAL(r8)    :: SBCX(imt,jmt,max_blocks_clinic),BBCX(imt,jmt,max_blocks_clinic)
      REAL(r8)    :: SBCY(imt,jmt,max_blocks_clinic),BBCY(imt,jmt,max_blocks_clinic)
!
      real(r8),allocatable,dimension(:,:) :: buffer

!     wjl 2021/08/03
!      real(r8),allocatable,dimension(:,:,:)::h0
!      real(r8),allocatable,dimension(:,:,:,:)::u,v
      real(r8),dimension(imt,jmt,max_blocks_clinic)::h0
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::u,v
!
!     ------------------------------------------------------------------
!     Pressure gradient
!     ------------------------------------------------------------------
     real(r8),dimension(:,:,:,:),allocatable::gg,dlu,dlv
     real(r8),dimension(:,:,:),allocatable::dlub,dlvb
      

end module dyn_mod

