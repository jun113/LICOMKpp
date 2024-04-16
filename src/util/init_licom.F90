SUBROUTINE init_licom

use precision_mod
use param_mod
use msg_mod

   IMPLICIT NONE

   call set_dyn_mod
   call set_forc_mod
   call set_pconst_mod
   call set_pmix_mod
   call allocate_work_mod
   call set_tracer_mod

END SUBROUTINE

SUBROUTINE allocate_work_mod
use work_mod
#include <../head/def-undef.h>
   IMPLICIT NONE
   integer :: error
   ! if (.not. allocated(pxb)) then
   !    allocate(pxb(imt, jmt, max_blocks_clinic), stat = ierr) 
   ! end if
   pxb = 0.0d0
   ! if (.not. allocated(pyb)) then
   !    allocate(pyb(imt, jmt, max_blocks_clinic), stat = ierr) 
   ! end if
   pyb = 0.0d0
   ! if (.not. allocated(pax)) then
   !    allocate(pax(imt, jmt, max_blocks_clinic), stat = ierr) 
   ! end if
   pax = 0.0d0
   ! if (.not. allocated(pay)) then
   !    allocate(pay(imt, jmt, max_blocks_clinic), stat = ierr) 
   ! end if
   pay = 0.0d0
   ! if (.not. allocated(whx)) then
   !    allocate(whx(imt, jmt, max_blocks_clinic), stat = ierr) 
   ! end if
   whx = 0.0d0
   ! if (.not. allocated(why)) then
   !    allocate(why(imt, jmt, max_blocks_clinic), stat = ierr) 
   ! end if
   why = 0.0d0
   ! if (.not. allocated(wgp)) then
   !    allocate(wgp(imt, jmt, max_blocks_clinic), stat = ierr) 
   ! end if
   wgp = 0.0d0

   ! if (.not. allocated(wkk)) then
   !    allocate(wkk(kmp1), stat = ierr) 
   ! end if
   wkk = 0.0d0

END SUBROUTINE allocate_work_mod

subroutine set_dyn_mod
use dyn_mod
#include <../head/def-undef.h>
   IMPLICIT NONE
   ub     = 0.0d0
   vb     = 0.0d0
   ubp    = 0.0d0
   vbp    = 0.0d0
   h0p    = 0.0d0
   up     = 0.0d0
   vp     = 0.0d0
   ws     = 0.0d0
   h0l    = 0.0d0
   h0f    = 0.0d0
   h0bl   = 0.0d0
   h0bf   = 0.0d0
   utl    = 0.0d0
   utf    = 0.0d0
   vtl    = 0.0d0
   vtf    = 0.0d0
   SBCX   = 0.0d0
   SBCY   = 0.0d0
   BBCX   = 0.0d0
   BBCY   = 0.0d0
!  buffer = 0.0d0
   h0     = 0.0d0
   u      = 0.0d0
   v      = 0.0d0
!  gg     = 0.0d0
!  dlu    = 0.0d0
!  dlv    = 0.0d0
!  dlub   = 0.0d0
!  dlvb   = 0.0d0
end SUBROUTINE
subroutine set_forc_mod
#include <../head/def-undef.h>
use forc_mod
   IMPLICIT NONE
!  su3        = 0.0d0
!  sv3        = 0.0d0
!  psa3       = 0.0d0
!  tsa3       = 0.0d0
!  qar3       = 0.0d0
!  uva3       = 0.0d0
!  swv3       = 0.0d0
!  cld3       = 0.0d0
!  sss3       = 0.0d0
!  sst3       = 0.0d0
!  nswv3      = 0.0d0
!  dqdt3      = 0.0d0
!  chloro3    = 0.0d0
!  wspd3      = 0.0d0
!  wspdu3     = 0.0d0
!  wspdv3     = 0.0d0
!  lwv3       = 0.0d0
!  seaice3    = 0.0d0
!  rain3      = 0.0d0
!  snow3      = 0.0d0
!  runoff3    = 0.0d0

   su         = 0.0d0
   sv         = 0.0d0
   psa        = 0.0d0
   tsa        = 0.0d0
   sss        = 0.0d0
   swv        = 0.0d0
   uva        = 0.0d0
   qar        = 0.0d0
   cld        = 0.0d0
   ddd        = 0.0d0
   qqq        = 0.0d0
   sst        = 0.0d0
   nswv       = 0.0d0
   dqdt       = 0.0d0
   chloro     = 0.0d0
   lwv        = 0.0d0
   seaice     = 0.0d0
   rain       = 0.0d0
   snow       = 0.0d0
   fresh      = 0.0d0
   runoff     = 0.0d0
   lthf       = 0.0d0
   sshf       = 0.0d0

   ustar      = 0.0d0
   buoytur    = 0.0d0
   buoysol    = 0.0d0

!  su3_io     = 0.0d0
!  sv3_io     = 0.0d0
!  psa3_io    = 0.0d0
!  tsa3_io    = 0.0d0
!  qar3_io    = 0.0d0
!  uva3_io    = 0.0d0
!  swv3_io    = 0.0d0
!  cld3_io    = 0.0d0
!  sss3_io    = 0.0d0
!  sst3_io    = 0.0d0
!  nswv3_io   = 0.0d0
!  dqdt3_io   = 0.0d0
!  chloro3_io = 0.0d0
!  wspdu3_io  = 0.0d0
!  wspdv3_io  = 0.0d0
!  lwv3_io    = 0.0d0
!  seaice3_io = 0.0d0
!  rain3_io   = 0.0d0
!  snow3_io   = 0.0d0
!  runoff3_io = 0.0d0

   restore    = 0.0d0
   tsf        = 0.0d0
   ssf        = 0.0d0
!   wave_dis   = 0.0d0
!   t10        = 0.0d0
!   u10        = 0.0d0
!   v10        = 0.0d0
!   slp        = 0.0d0
!   q10        = 0.0d0
!   swhf       = 0.0d0
!   lwhf       = 0.0d0
!   precr      = 0.0d0
!   precs      = 0.0d0
!   rf         = 0.0d0
!   si         = 0.0d0
!   buffer3d   = 0.0d0
!   w3d        = 0.0d0

end subroutine

subroutine set_pconst_mod
#include <../head/def-undef.h>
use pconst_mod
   IMPLICIT NONE

   j_global     = 0
   i_global     = 0
   vit          = 0.0d0
   viv          = 0.0d0
  !  ahv_back     = 0.0d0
  !  na           = 0
#if (defined NETCDF) || (defined ALL)
   lon          = 0.0d0
   lat          = 0.0d0
   lon_o        = 0.0d0
   lat_o        = 0.0d0
   ulon_o       = 0.0d0
   ulat_o       = 0.0d0
   lev          = 0.0d0
   lev1         = 0.0d0
#endif
   ! s_lon        = 0.0d0
   ! s_lat        = 0.0d0
   zkt          = 0.0d0
   dzp          = 0.0d0
   odzp         = 0.0d0
   odzt         = 0.0d0
   zkp          = 0.0d0
   ebea         = 0.0d0
   ebeb         = 0.0d0
   ebla         = 0.0d0
   eblb         = 0.0d0
   epea         = 0.0d0
   epeb         = 0.0d0
   epla         = 0.0d0
   eplb         = 0.0d0
   rrd1         = 0.0d0
   rrd2         = 0.0d0

   ohbt         = 0.0d0
   ohbu         = 0.0d0
   dzph         = 0.0d0
   hbx          = 0.0d0
   hby          = 0.0d0
   snlat        = 0.0d0
!  i_start      = 0
!  j_start      = 0
   to           = 0.0d0
   so           = 0.0d0
   c            = 0.0d0
   po           = 0.0d0
   akmu         = 0.0d0
   akmt         = 0.0d0
   akt          = 0.0d0
   am           = 0.0d0
   ah           = 0.0d0
   am3          = 0.0d0
  !  ah3          = 0.0d0
  !  amx          = 0.0d0
  !  amy          = 0.0d0
   nmonth       = 0
   nnmonth      = 0
#ifdef TIDEMIX
   fz_tide      = 0.0d0
   ak_tide      = 0.0d0
   fztidal      = 0.0d0
   richardson   = 0.0d0
   wp3_tidal    = 0.0d0
   ak_tide1     = 0.0d0
#endif
#ifdef CANUTOMIXOUT
   wp1_canuto   = 0.0d0
   wp2_canuto   = 0.0d0
   wp3_canuto   = 0.0d0
   wp4_canuto   = 0.0d0
   wp5_canuto   = 0.0d0
   wp6_canuto   = 0.0d0
   wp7_canuto   = 0.0d0
   wp8_canuto   = 0.0d0
   wp10_canuto  = 0.0d0
   wp11_canuto  = 0.0d0
   wp12_canuto  = 0.0d0
   wp13_canuto  = 0.0d0

   wk1_canuto   = 0.0d0
   wk2_canuto   = 0.0d0
   wk3_canuto   = 0.0d0
   wk4_canuto   = 0.0d0

   fcor_canuto  = 0.0d0
   fcort_canuto = 0.0d0

   alpha_canuto = 0.0d0
   beta_canuto  = 0.0d0
#endif
end SUBROUTINE

subroutine set_pmix_mod 
#include <../head/def-undef.h>
use pmix_mod
   IMPLICIT NONE

   rtst           = 0.0d0
   rtend          = 0.0d0
   rust           = 0.0d0
   ruend          = 0.0d0

   wndmix         = 0.0d0
   fricmx         = 0.0d0

   diff_cbt_back  = 0.0d0
   diff_cbt_limit = 0.0d0

   visc_cbu_back  = 0.0d0
   visc_cbu_limit = 0.0d0

!  ric            = 0.0d0
!  rict           = 0.0d0
!  rict_replace   = 0.0d0

!  rict_ref       = 0.0d0

!  rit            = 0.0d0
!  riu            = 0.0d0

   ricdt          = 0.0d0
   ricdttms       = 0.0d0

   ridt           = 0.0d0

   s2u            = 0.0d0
   s2t            = 0.0d0

#ifdef SOLAR
   pen            = 0.0d0
#endif

#ifdef SOLARCHLORO
   pen_chl        = 0.0d0
#endif
end subroutine

subroutine set_tracer_mod
#include <../head/def-undef.h>
use tracer_mod
   IMPLICIT NONE
   atb        = 0.0d0
   net        = 0.0d0

   at         = 0.0d0
!  restore_at = 0.0d0

   pdensity   = 0.0d0
   amld       = 0.0d0

   tend       = 0.0d0
   ax         = 0.0d0
   ay         = 0.0d0
   az         = 0.0d0

   dx         = 0.0d0
   dy         = 0.0d0
   dz         = 0.0d0

   penetrate  = 0.0d0

   dt_diff    = 0.0d0

   ddy        = 0.0d0

   dt_conv    = 0.0d0
   dt_restore = 0.0d0

#ifdef ISO
   aay_iso    = 0.0d0
   ddy_iso    = 0.0d0
   dx_iso     = 0.0d0

   dy_iso     = 0.0d0
#endif

   licomqice  = 0.0d0

   fw_norm2   = 0.0d0

end subroutine
