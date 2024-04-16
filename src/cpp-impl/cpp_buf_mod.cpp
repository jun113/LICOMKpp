#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/cpp_param_mod.h"
#include "../head/fortran_buf_mod.h"

namespace CppBufMod {

using CppParamMod::IMT;
using CppParamMod::JMT;

constexpr int NIBUFF = 100;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
double (&ibuffs)[NIBUFF] = buf_mod_mp_ibuffs_;
double (&ibuffr)[NIBUFF] = buf_mod_mp_ibuffr_;

int &nxg = buf_mod_mp_nxg_;
int &nyg = buf_mod_mp_nyg_;

int &nx = buf_mod_mp_nx_;
int &ny = buf_mod_mp_ny_;

int &nv = buf_mod_mp_nv_;

#ifdef USE_OCN_CARBON

#endif // USE_OCN_CARBON

double (*(&ifrac))[JMT][IMT] = buf_mod_mp_ifrac_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
double (&ibuffs)[NIBUFF] = __buf_mod_MOD_ibuffs;
double (&ibuffr)[NIBUFF] = __buf_mod_MOD_ibuffr;

int &nxg = __buf_mod_MOD_nxg;
int &nyg = __buf_mod_MOD_nyg;

int &nx = __buf_mod_MOD_nx;
int &ny = __buf_mod_MOD_ny;

int &nv = __buf_mod_MOD_nv;

#ifdef USE_OCN_CARBON

#endif // USE_OCN_CARBON

double (*(&ifrac))[JMT][IMT] = __buf_mod_MOD_ifrac;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

} // namespace CppBufMod
#endif // LICOM_ENABLE_FORTRAN
