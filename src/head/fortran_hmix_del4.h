#include "cpp_param_mod.h"
#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_HMIX_DEL4_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_HMIX_DEL4_H_
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppParamMod::MAX_BLOCKS_CLINIC;
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern double hmix_del4_mp_del4_dtn_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dts_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dte_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dtw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_duc_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dun_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dus_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_due_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_duw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dmc_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dmn_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dms_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dme_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dmw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_dum_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_ahf_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_amf_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del4_mp_del4_ratio_dxy_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern double hmix_del4_mp_ah_;
extern double hmix_del4_mp_am_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern double __hmix_del4_MOD_del4_dtn[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dts[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dte[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dtw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_duc[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dun[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dus[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_due[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_duw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dmc[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dmn[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dms[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dme[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dmw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_dum[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_ahf[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_amf[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del4_MOD_del4_ratio_dxy[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];

extern double __hmix_del4_MOD_ah;
extern double __hmix_del4_MOD_am;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_HMIX_DEL4_H_
