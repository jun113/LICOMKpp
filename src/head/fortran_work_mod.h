#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_WORK_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_WORK_MOD_H_

#include "cpp_param_mod.h"

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMP1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NTRA;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern double work_mod_mp_pxb_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double work_mod_mp_pyb_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double work_mod_mp_pax_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double work_mod_mp_pay_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double work_mod_mp_whx_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double work_mod_mp_why_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double work_mod_mp_wgp_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double work_mod_mp_wka_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

// extern double (*work_mod_mp_work_1_)[KM][JMT][IMT];
// extern double (*work_mod_mp_work_2_)[KM][JMT][IMT];
// extern double (*work_mod_mp_work_3_)[KM+1][JMT][IMT];

// extern double (*work_mod_mp_temp_)[KM][JMT][IMT];

// extern double (*work_mod_mp_uk_)[KM][JMT][IMT];
// extern double (*work_mod_mp_vk_)[KM][JMT][IMT];

extern double work_mod_mp_work_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double work_mod_mp_wkk_[KMP1];

// extern double (*work_mod_mp_wkb_)[KM][JMT][IMT];
// extern double (*work_mod_mp_wkc_)[KM][JMT][IMT];
// extern double (*work_mod_mp_wkd_)[KM][JMT][IMT];

// extern double (*work_mod_mp_tf_)[KM][JMT][IMT];
// extern double (*work_mod_mp_stf_)[JMT][IMT];

// extern float (*work_mod_mp_buffer_real4_)[IMT_GLOBAL];

// extern double (*work_mod_mp_work1_g_)[IMT_GLOBAL];
// extern double (*work_mod_mp_work2_g_)[IMT_GLOBAL];
// extern double (*work_mod_mp_work3_g_)[IMT_GLOBAL];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern double __work_mod_MOD_pxb[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __work_mod_MOD_pyb[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __work_mod_MOD_pax[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __work_mod_MOD_pay[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __work_mod_MOD_whx[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __work_mod_MOD_why[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __work_mod_MOD_wgp[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __work_mod_MOD_wka[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

// extern double (*__work_mod_MOD_work_1)[KM][JMT][IMT];
// extern double (*__work_mod_MOD_work_2)[KM][JMT][IMT];
// extern double (*__work_mod_MOD_work_3)[KM+1][JMT][IMT];

// extern double (*__work_mod_MOD_temp)[KM][JMT][IMT];

// extern double (*__work_mod_MOD_uk)[KM][JMT][IMT];
// extern double (*__work_mod_MOD_vk)[KM][JMT][IMT];

extern double __work_mod_MOD_work[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __work_mod_MOD_wkk[KMP1];

// extern double (*__work_mod_MOD_wkb)[KM][JMT][IMT];
// extern double (*__work_mod_MOD_wkc)[KM][JMT][IMT];
// extern double (*__work_mod_MOD_wkd)[KM][JMT][IMT];

// extern double (*__work_mod_MOD_tf)[KM][JMT][IMT];

// extern double (*__work_mod_MOD_stf)[JMT][IMT];

// extern float (*__work_mod_MOD_buffer_real4)[IMT_GLOBAL];

// extern double (*__work_mod_MOD_work1_g)[IMT_GLOBAL];
// extern double (*__work_mod_MOD_work2_g)[IMT_GLOBAL];
// extern double (*__work_mod_MOD_work3_g)[IMT_GLOBAL];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_WORK_MOD_H_