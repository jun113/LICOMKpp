#ifndef LICOM3_KOKKOSR_SRC_HEAD_FORTRAN_ISOPYC_MOD_H_
#define LICOM3_KOKKOSR_SRC_HEAD_FORTRAN_ISOPYC_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"

#ifdef ISO

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMP1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;

constexpr int NRPL = 5;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern double isopyc_mod_mp_dptlim_[NRPL+1];

extern double isopyc_mod_mp_fzisop_[KM];

extern double isopyc_mod_mp_slmxr_;

extern double isopyc_mod_mp_ahisop_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double isopyc_mod_mp_athkdf_[MAX_BLOCKS_CLINIC][JMT][IMT];

// extern double (*isopyc_mod_mp_e_)[3][JMT][KMP1][IMT];
// extern double (*isopyc_mod_mp_rhoi_)[NRPL][JMT][KM+1][IMT];

// extern double (*isopyc_mod_mp_k1_)[1][JMT][KM+1][IMT];
// extern double (*isopyc_mod_mp_k2_)[1][JMT][KM+1][IMT];
// extern double (*isopyc_mod_mp_k3_)[3][JMT][KM+1][IMT];

// extern double (*isopyc_mod_mp_adv_vetiso_)[JMT][KM][IMT];

extern double isopyc_mod_mp_adv_vbtiso_[MAX_BLOCKS_CLINIC][JMT][KM+1][IMT];
extern double isopyc_mod_mp_adv_vntiso_[MAX_BLOCKS_CLINIC][JMT][KM][IMT];


#ifdef isopycmixspatialvar
extern float isopyc_mod_mp_dciso1_[MAX_BLOCKS_CLINIC][JMT][KM][IMT];
extern float isopyc_mod_mp_dciso2_[MAX_BLOCKS_CLINIC][JMT][KM][IMT];

extern float isopyc_mod_mp_dslope_;
extern float isopyc_mod_mp_slopec_;
#endif // isopycmixspatialvar

extern double isopyc_mod_mp_kisrpl_[KM];

extern int isopyc_mod_mp_krplin_[NRPL];

extern double isopyc_mod_mp_zt_[KM];

extern double isopyc_mod_mp_dzw_[KM+1];
extern double isopyc_mod_mp_dzwr_[KM+1];

extern double isopyc_mod_mp_dzr_[KM];

extern double isopyc_mod_mp_tmask_[MAX_BLOCKS_CLINIC][JMT][KM][IMT];
extern double isopyc_mod_mp_f3_[MAX_BLOCKS_CLINIC][JMT][IMT];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern double __isopyc_mod_MOD_dptlim[NRPL+1];

extern double __isopyc_mod_MOD_fzisop[KM];

extern double __isopyc_mod_MOD_slmxr;

extern double __isopyc_mod_MOD_ahisop[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __isopyc_mod_MOD_athkdf[MAX_BLOCKS_CLINIC][JMT][IMT];

// extern double (*__isopyc_mod_MOD_e)[3][JMT][KMP1][IMT];
// extern double (*__isopyc_mod_MOD_rhoi)[NRPL][JMT][KM+1][IMT];

// extern double (*__isopyc_mod_MOD_k1)[1][JMT][KM+1][IMT];
// extern double (*__isopyc_mod_MOD_k2)[1][JMT][KM+1][IMT];
// extern double (*__isopyc_mod_MOD_k3)[3][JMT][KM+1][IMT];

// extern double (*__isopyc_mod_MOD_adv_vetiso)[JMT][KM][IMT];

extern double __isopyc_mod_MOD_adv_vbtiso[MAX_BLOCKS_CLINIC][JMT][KM+1][IMT];
extern double __isopyc_mod_MOD_adv_vntiso[MAX_BLOCKS_CLINIC][JMT][KM][IMT];


#ifdef isopycmixspatialvar
extern float __isopyc_mod_MOD_dciso1[MAX_BLOCKS_CLINIC][JMT][KM][IMT];
extern float __isopyc_mod_MOD_dciso2[MAX_BLOCKS_CLINIC][JMT][KM][IMT];

extern float __isopyc_mod_MOD_dslope;
extern float __isopyc_mod_MOD_slopec;
#endif // isopycmixspatialvar

extern double __isopyc_mod_MOD_kisrpl[KM];

extern int __isopyc_mod_MOD_krplin[NRPL];

extern double __isopyc_mod_MOD_zt[KM];

extern double __isopyc_mod_MOD_dzw[KM+1];
extern double __isopyc_mod_MOD_dzwr[KM+1];

extern double __isopyc_mod_MOD_dzr[KM];

extern double __isopyc_mod_MOD_tmask[MAX_BLOCKS_CLINIC][JMT][KM][IMT];
extern double __isopyc_mod_MOD_f3[MAX_BLOCKS_CLINIC][JMT][IMT];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
#endif // ISO
#endif // LICOM3_KOKKOSR_SRC_HEAD_FORTRAN_ISOPYC_MOD_H_
