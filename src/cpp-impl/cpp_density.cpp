#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_domain.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_tracer_mod.h"

void density () {
	using CppDomain::nblocks_clinic;
	using CppParamMod::KM;
	using CppParamMod::JST;
	using CppParamMod::JET;
	using CppParamMod::IMT;
  using CppPconstMod::c;
  using CppPconstMod::to;
  using CppPconstMod::so;
  using CppPconstMod::po;
  using CppPconstMod::vit;
  using CppTracerMod::at;
  using CppTracerMod::pdensity;

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = JST-1; j < JET; ++j) {
        for (int i = 0; i < IMT; ++i) {
          if (vit[iblock][k][j][i] > 0.0) {
            const double tq = at[iblock][0][k][j][i] - to[k]; 
            const double sq = at[iblock][1][k][j][i] - so[k]; 

            pdensity[iblock][k][j][i] = 1.0e+3 + po[k] +
	              (c[0][k] + (c[3][k] + c[6][k] * sq) *sq +
	              (c[2][k] +  c[7][k] * sq + c[5][k] * tq) * tq) * tq +
	              (c[1][k] + (c[4][k] + c[8][k] * sq) * sq) * sq;
          } else {
            pdensity[iblock][k][j][i] = 0.0;
          }
        }
      }
    }
  }
	return ;
}
#endif // LICOM_ENABLE_FORTRAN