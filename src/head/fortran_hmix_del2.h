#include "cpp_param_mod.h"
#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_HMIX_DEL2_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_HMIX_DEL2_H_
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
using CppParamMod::MAX_BLOCKS_CLINIC;
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern double hmix_del2_mp_dtn_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dts_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dte_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dtw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_duc_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dun_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dus_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_due_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_duw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dmc_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dmn_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dms_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dme_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dmw_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_dum_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_ahf_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_amf_[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_ah_;
extern double hmix_del2_mp_am_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern double __hmix_del2_MOD_dtn[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dts[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dte[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dtw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_duc[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dun[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dus[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_due[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_duw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dmc[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dmn[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dms[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dme[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dmw[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_dum[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_ahf[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_amf[MAX_BLOCKS_CLINIC][NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_ah;
extern double __hmix_del2_MOD_am;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_HMIX_DEL2_H_

