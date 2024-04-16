#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_pconst_mod.h"
double dens (const double &tq, const double &sq, const int &kk) {
  double dens;
  using CppPconstMod::c;
  dens = (c[0][kk] + (c[3][kk] + c[6][kk] * sq) * sq +
         (c[2][kk] +  c[7][kk] * sq + c[5][kk] * tq) * tq) * tq +
         (c[1][kk] + (c[4][kk] + c[8][kk] * sq) * sq) * sq;
  return dens;
}
#endif // LICOM_ENABLE_FORTRAN