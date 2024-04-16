#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_TRACER_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_TRACER_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NTRA;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern double tracer_mod_mp_atb_[MAX_BLOCKS_CLINIC][NTRA][KM+1][JMT][IMT];
extern double tracer_mod_mp_net_[MAX_BLOCKS_CLINIC][NTRA][JMT][IMT];

extern double tracer_mod_mp_at_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
// extern double (*tracer_mod_mp_restore_at_)[NTRA][JMT][IMT];

extern double tracer_mod_mp_pdensity_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double tracer_mod_mp_amld_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double tracer_mod_mp_tend_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_ax_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_ay_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_az_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double tracer_mod_mp_dx_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_dy_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_dz_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double tracer_mod_mp_penetrate_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double tracer_mod_mp_dt_diff_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double tracer_mod_mp_ddy_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double tracer_mod_mp_dt_conv_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_dt_restore_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

#ifdef ISO
extern double tracer_mod_mp_aay_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_ddy_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double tracer_mod_mp_ax_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_ay_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_az_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double tracer_mod_mp_dx_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_dy_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double tracer_mod_mp_dz_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
#endif // ISO

extern double tracer_mod_mp_licomqice_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double tracer_mod_mp_fw_norm2_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern double __tracer_mod_MOD_atb[MAX_BLOCKS_CLINIC][NTRA][KM+1][JMT][IMT];
extern double __tracer_mod_MOD_net[MAX_BLOCKS_CLINIC][NTRA][JMT][IMT];

extern double __tracer_mod_MOD_at[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
// extern double (*__tracer_mod_MOD_restore_at)[NTRA][JMT][IMT];

extern double __tracer_mod_MOD_pdensity[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __tracer_mod_MOD_amld[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __tracer_mod_MOD_tend[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_ax[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_ay[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_az[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double __tracer_mod_MOD_dx[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_dy[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_dz[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double __tracer_mod_MOD_penetrate[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double __tracer_mod_MOD_dt_diff[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double __tracer_mod_MOD_ddy[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double __tracer_mod_MOD_dt_conv[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_dt_restore[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

#ifdef ISO
extern double __tracer_mod_MOD_aay_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_ddy_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double __tracer_mod_MOD_ax_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_ay_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_az_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double __tracer_mod_MOD_dx_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_dy_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __tracer_mod_MOD_dz_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
#endif // ISO

extern double __tracer_mod_MOD_licomqice[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __tracer_mod_MOD_fw_norm2;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_TRACER_MOD_H_
