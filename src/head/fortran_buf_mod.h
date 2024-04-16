#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_BUF_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_BUF_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"

using CppParamMod::IMT;
using CppParamMod::JMT;

constexpr int NIBUFF = 100;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern double buf_mod_mp_ibuffs_[NIBUFF];
extern double buf_mod_mp_ibuffr_[NIBUFF];

extern int buf_mod_mp_nxg_;
extern int buf_mod_mp_nyg_;

extern int buf_mod_mp_nx_;
extern int buf_mod_mp_ny_;

extern int buf_mod_mp_nv_;



#ifdef USE_OCN_CARBON

#endif

extern double (*buf_mod_mp_ifrac_)[JMT][IMT];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern double __buf_mod_MOD_ibuffs[NIBUFF];
extern double __buf_mod_MOD_ibuffr[NIBUFF];

extern int __buf_mod_MOD_nxg;
extern int __buf_mod_MOD_nyg;

extern int __buf_mod_MOD_nx;
extern int __buf_mod_MOD_ny;

extern int __buf_mod_MOD_nv;



#ifdef USE_OCN_CARBON

#endif

extern double (*__buf_mod_MOD_ifrac)[JMT][IMT];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_BUF_MOD_H_
