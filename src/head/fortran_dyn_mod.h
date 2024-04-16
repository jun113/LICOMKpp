#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_DYN_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_DYN_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::KMP1;
using CppParamMod::IMT_GLOBAL;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern double dyn_mod_mp_ub_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_vb_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_ubp_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_vbp_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_h0p_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double dyn_mod_mp_up_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double dyn_mod_mp_vp_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double dyn_mod_mp_ws_[MAX_BLOCKS_CLINIC][KMP1][JMT][IMT];

extern double dyn_mod_mp_h0l_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_h0f_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_h0bl_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_h0bf_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double dyn_mod_mp_utl_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double dyn_mod_mp_utf_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double dyn_mod_mp_vtl_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double dyn_mod_mp_vtf_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double dyn_mod_mp_sbcx_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_bbcx_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_sbcy_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double dyn_mod_mp_bbcy_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (*dyn_mod_mp_buffer_)[IMT_GLOBAL];

extern double dyn_mod_mp_h0_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double dyn_mod_mp_u_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double dyn_mod_mp_v_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

// extern double (*dyn_mod_mp_gg_)[KM][JMT][IMT];
// extern double (*dyn_mod_mp_dlu_)[KM][JMT][IMT];
// extern double (*dyn_mod_mp_dlv_)[KM][JMT][IMT];
// extern double (*dyn_mod_mp_dlub_)[JMT][IMT];
// extern double (*dyn_mod_mp_dlvb_)[JMT][IMT];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern double __dyn_mod_MOD_ub[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_vb[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_ubp[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_vbp[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_h0p[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __dyn_mod_MOD_up[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __dyn_mod_MOD_vp[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double __dyn_mod_MOD_ws[MAX_BLOCKS_CLINIC][KMP1][JMT][IMT];

extern double __dyn_mod_MOD_h0l[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_h0f[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_h0bl[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_h0bf[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __dyn_mod_MOD_utl[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __dyn_mod_MOD_utf[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __dyn_mod_MOD_vtl[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __dyn_mod_MOD_vtf[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double __dyn_mod_MOD_sbcx[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_bbcx[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_sbcy[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __dyn_mod_MOD_bbcy[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double (*__dyn_mod_MOD_buffer)[IMT_GLOBAL];

extern double __dyn_mod_MOD_h0[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __dyn_mod_MOD_u[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __dyn_mod_MOD_v[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

// extern double (*__dyn_mod_MOD_gg)[KM][JMT][IMT];
// extern double (*__dyn_mod_MOD_dlu)[KM][JMT][IMT];
// extern double (*__dyn_mod_MOD_dlv)[KM][JMT][IMT];
// extern double (*__dyn_mod_MOD_dlub)[JMT][IMT];
// extern double (*__dyn_mod_MOD_dlvb)[JMT][IMT];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_DYN_MOD_H_
