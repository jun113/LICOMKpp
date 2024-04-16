#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_domain.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_tracer_mod.h"
#include "../head/cpp_shr_const_mod.h"

#include <cmath>

void cpp_icesnow(){

  using CppGrid::kmt;
  using CppParamMod::KM;  
  using CppParamMod::JMT;
  using CppParamMod::IMT;
  using CppPconstMod::cp;
  using CppPconstMod::tbice;
  using CppPconstMod::dzp;
  using CppTracerMod::at;
  using CppTracerMod::licomqice;
  using CppDomain::nblocks_clinic;
  using CppShrConstMod::SHR_CONST_LATICE; 

  const double sal_ice = 4.0;
  const double sal_ocn = 34.7;
  const double heat_ice_fusion = SHR_CONST_LATICE;

  double tdiff;
   
 //icesnow 1
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0 ; k < 1; ++k) {   // useless, just for consistency
      for (int j = 0 ; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          if(at[iblock][0][k][j][i] < tbice && kmt[iblock][j][i] > 0) {
            tdiff = tbice - at[iblock][0][k][j][i];
            // useless, just for consistency
            licomqice[iblock][j][i] += tdiff * dzp[k] / dzp[0];
            at[iblock][1][k][j][i]  += tdiff * (sal_ocn - sal_ice) 
                * cp / heat_ice_fusion * 0.001;
            at[iblock][0][k][j][i] = tbice;
          }//end if
        }
      }
    }
    for (int j = 0 ; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        if(licomqice[iblock][j][i] > 0.0 && 
            at[iblock][0][0][j][i] > tbice) {

          tdiff = fmin((at[iblock][0][0][j][i] - tbice), 
              licomqice[iblock][j][i]);

          licomqice[iblock][j][i] -= tdiff;
          at[iblock][0][0][j][i]  -= tdiff;
          at[iblock][1][0][j][i]  -= tdiff * (sal_ocn - sal_ice) 
              * cp / heat_ice_fusion * 0.001;
        }//end if
      }
    }
    for (int k = 0 ; k < KM; ++k) {
      for (int j = 0 ; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          if (at[iblock][0][k][j][i] < tbice) {
            at[iblock][0][k][j][i] = tbice; 
          }  
        }
      }
    }//end for
  }//end icesnow 1
  return ;
} // end cpp_icesnow
#endif // LICOM_ENABLE_FORTRAN
