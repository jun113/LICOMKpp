#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_domain.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"

using CppDomain::nblocks_clinic;
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;

void vinteg(const double (&wk3)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
    double (&wk2)[MAX_BLOCKS_CLINIC][JMT][IMT]) {

  using CppPconstMod::dzp;
  using CppPconstMod::viv;
  using CppPconstMod::ohbu;

  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wk2[iblock][j][i] = 0.0;
      }
    }
  }
  for (int iblock = 0; iblock < nblocks_clinic; ++iblock) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 0; i < IMT; ++i) {
          wk2[iblock][j][i] = wk2[iblock][j][i] + dzp[k] *
              ohbu[iblock][j][i] * wk3[iblock][k][j][i] *
              viv[iblock][k][j][i];
        }
      }
    }
  }
  return ;
}

#endif // LICOM_ENABLE_FORTRAN