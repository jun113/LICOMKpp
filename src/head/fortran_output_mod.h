#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_OUTPUT_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_OUTPUT_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NTRA;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
#ifdef DAILYACC
extern float output_mod_mp_z0daily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_mlddaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_lthfdaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_sshfdaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_lwvdaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_swvdaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_precdaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_evapdaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_roffdaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_sudaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_svdaily_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_tsdaily_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float output_mod_mp_ssdaily_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float output_mod_mp_usdaily_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float output_mod_mp_vsdaily_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float output_mod_mp_wsdaily_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#ifdef LOWRES
extern float output_mod_mp_z0mon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_himon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_hdmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_qicemon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_lthfmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_sshfmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_lwvmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_swvmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_sumon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_svmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_runoffmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_freshmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_wsmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_tsmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_ssmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_usmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_vsmon_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float output_mod_mp_icmon_[MAX_BLOCKS_CLINIC][2][JMT][IMT];

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
extern float output_mod_mp_athkdfmon_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // ISO_TYPE_BF

#ifdef ISO
extern float output_mod_mp_axmon_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float output_mod_mp_aymon_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float output_mod_mp_azmon_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float output_mod_mp_dxmon_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float output_mod_mp_dymon_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float output_mod_mp_dzmon_iso_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float output_mod_mp_wsmon_iso_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float output_mod_mp_vsmon_iso_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#ifdef ISOOUT
extern float output_mod_mp_azmon_vetisomon_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float output_mod_mp_azmon_vntisomon_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float output_mod_mp_azmon_vbtisomon_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
extern float output_mod_mp_azmon_a3mon_[KM][JMT][IMT];
#endif // SMAG_OUT

extern double output_mod_mp_err_norm2mon_;
#endif // LOWRES
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
#ifdef DAILYACC
extern float __output_mod_MOD_z0daily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_mlddaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_lthfdaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_sshfdaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_lwvdaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_swvdaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_precdaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_evapdaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_roffdaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_sudaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_svdaily[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_tsdaily[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float __output_mod_MOD_ssdaily[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float __output_mod_MOD_usdaily[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float __output_mod_MOD_vsdaily[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float __output_mod_MOD_wsdaily[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#ifdef DAILYBUGDET

#endif // DAILYBUGDET
#endif // DAILYACC

#ifdef LOWRES
extern float __output_mod_MOD_z0mon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_himon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_hdmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_qicemon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_lthfmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_sshfmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_lwvmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_swvmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_sumon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_svmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_runoffmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_freshmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_wsmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_tsmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_ssmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_usmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_vsmon[MAX_BLOCKS_CLINIC][JMT][IMT];
extern float __output_mod_MOD_icmon[MAX_BLOCKS_CLINIC][2][JMT][IMT];

#ifdef TIDEMIX

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

#endif // CANUTOMIXOUT

#ifdef ISO_TYPE_BF
extern float __output_mod_MOD_athkdfmon[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // ISO_TYPE_BF

#ifdef ISO
extern float __output_mod_MOD_axmon_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float __output_mod_MOD_aymon_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float __output_mod_MOD_azmon_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float __output_mod_MOD_dxmon_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float __output_mod_MOD_dymon_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float __output_mod_MOD_dzmon_iso[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern float __output_mod_MOD_wsmon_iso[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float __output_mod_MOD_vsmon_iso[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#ifdef ISOOUT
extern float __output_mod_MOD_azmon_vetisomon[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float __output_mod_MOD_azmon_vntisomon[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern float __output_mod_MOD_azmon_vbtisomon[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // ISOOUT
#endif // ISO

#ifdef SMAG_OUT
extern float __output_mod_MOD_azmon_a3mon[KM][JMT][IMT];
#endif // SMAG_OUT

extern double __output_mod_MOD_err_norm2mon;
#endif // LOWRES
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_OUTPUT_MOD_H_
