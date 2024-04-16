#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_BUF_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_BUF_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"

namespace CppBufMod {

using CppParamMod::IMT;
using CppParamMod::JMT;

constexpr int NIBUFF = 100;
constexpr int NSND   = 7;
constexpr int NRCV   = 14;

extern double (&ibuffs)[NIBUFF];
extern double (&ibuffr)[NIBUFF];

extern int &nxg;
extern int &nyg;

extern int &nx;
extern int &ny;

extern int &nv;

#ifdef USE_OCN_CARBON

#endif // USE_OCN_CARBON

extern double* (&ifrac)[JMT][IMT];

} // namespace CppBufMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_BUF_MOD_H_
