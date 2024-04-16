#include "cpp_param_mod.h"
#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_HMIX_DEL2_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_HMIX_DEL2_H_
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern double (*hmix_del2_mp_dtn_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dts_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dte_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dtw_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_duc_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dun_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dus_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_due_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_duw_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dmc_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dmn_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dms_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dme_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dmw_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_dum_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_ahf_)[NY_BLOCK][NX_BLOCK];
extern double (*hmix_del2_mp_amf_)[NY_BLOCK][NX_BLOCK];
extern double hmix_del2_mp_ah_;
extern double hmix_del2_mp_am_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern double (*__hmix_del2_MOD_dtn)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dts)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dte)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dtw)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_duc)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dun)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dus)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_due)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_duw)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dmc)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dmn)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dms)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dme)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dmw)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_dum)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_ahf)[NY_BLOCK][NX_BLOCK];
extern double (*__hmix_del2_MOD_amf)[NY_BLOCK][NX_BLOCK];
extern double __hmix_del2_MOD_ah;
extern double __hmix_del2_MOD_am;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_HMIX_DEL2_H_

