#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_domain.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"

using CppParamMod::KM;  
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::MAX_BLOCKS_CLINIC;

void invtrit(double (&wk)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
    const double (&topbc)[MAX_BLOCKS_CLINIC][JMT][IMT], 
				const double (&dcb)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
        		const double &aidif, const double &c2dtts) {
    
	using CppDomain::nblocks_clinic;
  using CppGrid::kmt;
  using CppPconstMod::vit;
  using CppPconstMod::odzp;
  using CppPconstMod::odzt;

  double a8[KM], b8[KM], c8[KM], d8[KM];
  double e8[KM+1], f8[KM+1];

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 2; j < JMT-2; ++j) {
      for (int i = 2; i < IMT-2; ++i) {

        if (kmt[iblock][j][i] > 0) {
          const int kz = kmt[iblock][j][i] - 1;

          for (int k = 1; k <= kz; ++k) {
            a8[k] = dcb[iblock][k-1][j][i] * odzt[k  ] * odzp[k] * c2dtts * aidif;
            d8[k] = wk[iblock][k][j][i];
          }

          for (int k = 1; k <= kz-1; ++k) {
            c8[k] = dcb[iblock][k  ][j][i] * odzt[k+1] * odzp[k] * c2dtts * aidif;
            b8[k] = 1.0 + a8[k] + c8[k];
            e8[k] = 0.0;
            f8[k] = 0.0;
          }
          // k = 0
          a8[0] = odzp[0] * c2dtts * aidif;
          c8[0] = dcb[iblock][0][j][i] * odzt[1] * odzp[0] * c2dtts * aidif;
          b8[0] = 1.0 + c8[0];
          d8[0] = wk[iblock][0][j][i];
          e8[0] = 0.0;
          f8[0] = 0.0;

          b8[kz] = 1.0 + a8[kz];
          c8[kz] = odzp[kz] * c2dtts * aidif;

          e8[kz+1] = 0.0;
          f8[kz+1] = 0.0;

          for (int k = kz; k >= 0; --k) {
            const double g0 = 1.0 / (b8[k] - c8[k] * e8[k+1]);
            e8[k] = a8[k] * g0;
            f8[k] = (d8[k] + c8[k] * f8[k+1]) * g0;
          }

          wk[iblock][0][j][i] = (e8[0] * topbc[iblock][j][i] + f8[0])
              * vit[iblock][0][j][i];

          for (int k = 1; k <= kz; ++k) {
            wk[iblock][k][j][i] = (e8[k] * wk[iblock][k-1][j][i]
                + f8[k]) * vit[iblock][k][j][i];
          }
        }
      }
    }
  }
  return ;
}

void invtriu(
    double (&wk)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
    const double (&topbc)[MAX_BLOCKS_CLINIC][JMT][IMT],
    const double (&bombc)[MAX_BLOCKS_CLINIC][JMT][IMT],
    const double (&dcb)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
    const double &aidif,
    const double &c2dtc){

	using CppDomain::nblocks_clinic;
	using CppGrid::kmu;//M NY NX
	using CppPconstMod::odzp;//K
	using CppPconstMod::odzt;//K
	using CppPconstMod::viv;//M K J I

  double a8[KM], b8[KM], c8[KM], d8[KM], e8[KM+1], f8[KM+1];
  int k, kz;
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 2; j < JMT - 2; ++j) {
      for (int i = 2; i < IMT - 2; ++i) {
        if (kmu[iblock][j][i] > 0) {
          kz = kmu[iblock][j][i]; 
          for (k = 1; k < kz ; ++k) { 
            a8[k] = dcb[iblock][k-1][j][i] 
                * odzt[k] * odzp[k] 
                    * c2dtc * aidif;
            d8[k] = wk[iblock][k][j][i];
          }
          for (k = 1; k < kz - 1; ++k) { 
            c8[k] = dcb[iblock][k][j][i] 
                * odzt[k+1] * odzp[k] 
                    * c2dtc * aidif;
            b8[k] = 1.0 +  a8[k] + c8[k];

            e8[k] = 0.0;
            f8[k] = 0.0;
          }
          //B.C. AT TOP
          k = 0;
          a8[k] = odzp[k] * c2dtc * aidif;
          c8[k] = dcb[iblock][k][j][i] 
              * odzt[k+1] * odzp[k] 
                  * c2dtc * aidif;
          b8[k] = 1.0 + c8[k];
          d8[k] = wk[iblock][k][j][i];
          e8[k] = 0.0;
          f8[k] = 0.0;
          //B.C. AT BOTTOM
          b8[kz-1] = 1.0 + a8[kz-1];
          c8[kz-1] = odzp[kz-1] * c2dtc * aidif;
          e8[kz] = 0.0;
          f8[kz] = 0.0;
          d8[kz-1] = wk[iblock][kz-1][j][i] 
              - bombc[iblock][j][i] 
                  * odzp[kz-1] * c2dtc * aidif;
          //NOW INVERT
          for (k = kz - 1; k >= 0; --k) {
            double g0 = 1.0 / (b8[k] - c8[k] * e8[k+1]);
            e8[k] = a8[k] * g0;
            f8[k] = (d8[k] + c8[k] * 
                f8[k+1]) * g0;
          }
         //B.C. AT SURFACE 
          wk[iblock][0][j][i] = (e8[0] * topbc[iblock][j][i] + f8[0]) 
              * viv[iblock][0][j][i];
          for (k = 1; k < kz ; ++k) {
            wk[iblock][k][j][i] = (e8[k] * wk[iblock][k-1][j][i] + f8[k]) 
                * viv[iblock][k][j][i];
          }
        }//end if
      }
    }
  }//end for
  return ;
}//end invtru  
#endif // LICOM_ENABLE_FORTRAN