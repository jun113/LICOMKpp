#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_PMIX_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_PMIX_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMM1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL

extern int pmix_mod_mp_rtst_;
extern int pmix_mod_mp_rtend_;
extern int pmix_mod_mp_rust_;
extern int pmix_mod_mp_ruend_;

extern double pmix_mod_mp_wndmix_;
extern double pmix_mod_mp_fricmx_;

extern double pmix_mod_mp_diff_cbt_back_;
extern double pmix_mod_mp_diff_cbt_limit_;

extern double pmix_mod_mp_visc_cbu_back_;
extern double pmix_mod_mp_visc_cbu_limit_;

// extern double (*pmix_mod_mp_ric_)[KMM1][JMT][IMT];
// extern double (*pmix_mod_mp_rict_)[KMM1][JMT][IMT];
// extern double (*pmix_mod_mp_rict_replace_)[KMM1][JMT][IMT];

// extern double (*pmix_mod_mp_rict_ref_)[JMT][IMT];

// extern double (*pmix_mod_mp_rit_)[KMM1][JMT][IMT];
// extern double (*pmix_mod_mp_riu_)[KM+1][JMT][IMT];

extern double pmix_mod_mp_ricdt_[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
extern double pmix_mod_mp_ricdttms_[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

extern double pmix_mod_mp_ridt_[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

extern double pmix_mod_mp_ridt_[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

extern double pmix_mod_mp_s2u_[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
extern double pmix_mod_mp_s2t_[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

#ifdef SOLAR
extern double pmix_mod_mp_pen_[KMM1];
#endif // SOLAR

#ifdef SOLARCHLORO
extern double pmix_mod_mp_pen_chl_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // SOLARCHLORO
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern int __pmix_mod_MOD_rtst;
extern int __pmix_mod_MOD_rtend;
extern int __pmix_mod_MOD_rust;
extern int __pmix_mod_MOD_ruend;

extern double __pmix_mod_MOD_wndmix;
extern double __pmix_mod_MOD_fricmx;

extern double __pmix_mod_MOD_diff_cbt_back;
extern double __pmix_mod_MOD_diff_cbt_limit;

extern double __pmix_mod_MOD_visc_cbu_back;
extern double __pmix_mod_MOD_visc_cbu_limit;

// extern double (*__pmix_mod_MOD_ric)[KMM1][JMT][IMT];
// extern double (*__pmix_mod_MOD_rict)[KMM1][JMT][IMT];
// extern double (*__pmix_mod_MOD_rict_replace)[KMM1][JMT][IMT];

// extern double (*__pmix_mod_MOD_rict_ref)[JMT][IMT];

// extern double (*__pmix_mod_MOD_rit)[KMM1][JMT][IMT];
// extern double (*__pmix_mod_MOD_riu)[KM+1][JMT][IMT];

extern double __pmix_mod_MOD_ricdt[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
extern double __pmix_mod_MOD_ricdttms[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

extern double __pmix_mod_MOD_ridt[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

extern double __pmix_mod_MOD_ridt[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

extern double __pmix_mod_MOD_s2u[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];
extern double __pmix_mod_MOD_s2t[MAX_BLOCKS_CLINIC][KMM1][JMT][IMT];

#ifdef SOLAR
extern double __pmix_mod_MOD_pen[KMM1];
#endif // SOLAR

#ifdef SOLARCHLORO
extern double __pmix_mod_MOD_pen_chl[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // SOLARCHLORO
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_PMIX_MOD_H_
