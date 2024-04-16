SUBROUTINE pop_haloupdate_smts(errorCode)
use precision_mod, only: r8
use dyn_mod, only: vtl
use domain, only: POP_haloClinic
use POP_HaloMod, only: POP_HaloUpdate
use POP_GridHorzMod, only: POP_gridHorzLocCenter, POP_fieldKindScalar

!LPF20160515
implicit none
integer::errorCode

   call POP_HaloUpdate(VTL , POP_haloClinic, POP_gridHorzLocCenter,&
      POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

END SUBROUTINE
