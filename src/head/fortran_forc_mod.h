#ifndef LICOM3_KOKKOS_SRC_FORTRAN_FORC_MOD_H_
#define LICOM3_KOKKOS_SRC_FORTRAN_FORC_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NTRA;
using CppParamMod::S_IMT;
using CppParamMod::S_JMT;

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
// extern double (*forc_mod_mp_su3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_sv3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_psa3_)[12][JMT][IMT];
// extern double forc_mod_mp_tsa3_[MAX_BLOCKS_CLINIC][JMT][IMT];
// extern double (*forc_mod_mp_qar3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_uva3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_swv3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_cld3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_sss3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_sst3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_nswv3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_dqdt3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_chloro3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_wspd3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_wspdu3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_wspdv3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_lwv3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_seaice3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_rain3_)[12][JMT][IMT];
// extern double (*forc_mod_mp_snow3_)[12][JMT][IMT];

extern double forc_mod_mp_su_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_sv_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_psa_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_tsa_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_sss_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_swv_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_uva_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_qar_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_cld_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_ddd_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_qqq_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_sst_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_nswv_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_dqdt_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_chloro_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_lwv_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_seaice_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_rain_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_snow_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_fresh_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_runoff_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_lthf_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_sshf_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double forc_mod_mp_ustar_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_buoytur_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_buoysol_[MAX_BLOCKS_CLINIC][JMT][IMT];

/*
extern double *(forc_mod_mp_su3_io_)[];
extern double *(forc_mod_mp_sv3_io_)[];
extern double *(forc_mod_mp_psa3_io_)[];
extern double *(forc_mod_mp_tsa3_io_)[];
extern double *(forc_mod_mp_qar3_io_)[];
extern double *(forc_mod_mp_uva3_io_)[];
extern double *(forc_mod_mp_swv3_io_)[];
extern double *(forc_mod_mp_cld3_io_)[];
extern double *(forc_mod_mp_sss3_io_)[];
extern double *(forc_mod_mp_sst3_io_)[];
extern double *(forc_mod_mp_nswv3_io_)[];
extern double *(forc_mod_mp_dqdt3_io_)[];
extern double *(forc_mod_mp_chloro3_io_)[];
extern double *(forc_mod_mp_wspd3_io_)[];
extern double *(forc_mod_mp_wspdu3_io_)[];
extern double *(forc_mod_mp_wspdv3_io_)[];
extern double *(forc_mod_mp_lwv3_io_)[];
extern double *(forc_mod_mp_seaice_io_)[];
extern double *(forc_mod_mp_rain_io_)[];
extern double *(forc_mod_mp_snow_io_)[];
extern double *(forc_mod_mp_fresh_io_)[];
extern double *(forc_mod_mp_runoff_io_)[];
*/

extern double forc_mod_mp_restore_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double forc_mod_mp_tsf_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double forc_mod_mp_ssf_[MAX_BLOCKS_CLINIC][JMT][IMT];

#ifdef TIDEMIX
extern double forc_mod_mp_wave_dis_[MAX_BLOCKS_CLINIC][JMT][IMT];
#endif // TIDEMIX

// extern double (*forc_mod_mp_t10_)[S_JMT][S_IMT];
// extern double (*forc_mod_mp_u10_)[S_JMT][S_IMT];
// extern double (*forc_mod_mp_v10_)[S_JMT][S_IMT];
// extern double (*forc_mod_mp_slp_)[S_JMT][S_IMT];
// extern double (*forc_mod_mp_q10_)[S_JMT][S_IMT];
// extern double (*forc_mod_mp_swhf_)[S_JMT][S_IMT];
// extern double (*forc_mod_mp_lwhf_)[S_JMT][S_IMT];

// extern double (*forc_mod_mp_precr_)[S_JMT][S_IMT];
// extern double (*forc_mod_mp_precs_)[S_JMT][S_IMT];

// extern double (*forc_mod_mp_rf_)[S_JMT][S_IMT];
// extern double (*forc_mod_mp_si_)[S_JMT][S_IMT];
// extern double (*forc_mod_mp_buffer3d_)[KM][JMT][IMT];
// extern double (*forc_mod_mp_w3d_)[JMT_GLOBAL][IMT_GLOBAL];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
// extern double (*__forc_mod_MOD_su3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_sv3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_psa3)[12][JMT][IMT];
// extern double __forc_mod_MOD_tsa3[MAX_BLOCKS_CLINIC][JMT][IMT];
// extern double (*__forc_mod_MOD_qar3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_uva3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_swv3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_cld3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_sss3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_sst3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_nswv3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_dqdt3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_chloro3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_wspd3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_wspdu3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_wspdv3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_lwv3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_seaice3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_rain3)[12][JMT][IMT];
// extern double (*__forc_mod_MOD_snow3)[12][JMT][IMT];

extern double __forc_mod_MOD_su[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_sv[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_psa[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_tsa[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_sss[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_swv[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_uva[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_qar[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_cld[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_ddd[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_qqq[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_sst[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_nswv[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_dqdt[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_chloro[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_lwv[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_seaice[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_rain[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_snow[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_fresh[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_runoff[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_lthf[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_sshf[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __forc_mod_MOD_ustar[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_buoytur[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_buoysol[MAX_BLOCKS_CLINIC][JMT][IMT];

/*
extern double *(__forc_mod_MOD_su3_io)[];
extern double *(__forc_mod_MOD_sv3_io)[];
extern double *(__forc_mod_MOD_psa3_io)[];
extern double *(__forc_mod_MOD_tsa3_io)[];
extern double *(__forc_mod_MOD_qar3_io)[];
extern double *(__forc_mod_MOD_uva3_io)[];
extern double *(__forc_mod_MOD_swv3_io)[];
extern double *(__forc_mod_MOD_cld3_io)[];
extern double *(__forc_mod_MOD_sss3_io)[];
extern double *(__forc_mod_MOD_sst3_io)[];
extern double *(__forc_mod_MOD_nswv3_io)[];
extern double *(__forc_mod_MOD_dqdt3_io)[];
extern double *(__forc_mod_MOD_chloro3_io)[];
extern double *(__forc_mod_MOD_wspd3_io)[];
extern double *(__forc_mod_MOD_wspdu3_io)[];
extern double *(__forc_mod_MOD_wspdv3_io)[];
extern double *(__forc_mod_MOD_lwv3_io)[];
extern double *(__forc_mod_MOD_seaice_io)[];
extern double *(__forc_mod_MOD_rain_io)[];
extern double *(__forc_mod_MOD_snow_io)[];
extern double *(__forc_mod_MOD_fresh_io)[];
extern double *(__forc_mod_MOD_runoff_io)[];
*/

extern double __forc_mod_MOD_restore[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];
extern double __forc_mod_MOD_tsf[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __forc_mod_MOD_ssf[MAX_BLOCKS_CLINIC][JMT][IMT];

#ifdef TIDEMIX
extern double __forc_mod_MOD_wave_dis[MAX_BLOCKS_CLINIC][JMT][IMT];
#endif // TIDEMIX

// extern double (*__forc_mod_MOD_t10)[S_JMT][S_IMT];
// extern double (*__forc_mod_MOD_u10)[S_JMT][S_IMT];
// extern double (*__forc_mod_MOD_v10)[S_JMT][S_IMT];
// extern double (*__forc_mod_MOD_slp)[S_JMT][S_IMT];
// extern double (*__forc_mod_MOD_q10)[S_JMT][S_IMT];
// extern double (*__forc_mod_MOD_swhf)[S_JMT][S_IMT];
// extern double (*__forc_mod_MOD_lwhf)[S_JMT][S_IMT];

// extern double (*__forc_mod_MOD_precr)[S_JMT][S_IMT];
// extern double (*__forc_mod_MOD_precs)[S_JMT][S_IMT];

// extern double (*__forc_mod_MOD_rf)[S_JMT][S_IMT];
// extern double (*__forc_mod_MOD_si)[S_JMT][S_IMT];
// extern double (*__forc_mod_MOD_buffer3d)[KM][JMT][IMT];
// extern double (*__forc_mod_MOD_w3d)[JMT_GLOBAL][IMT_GLOBAL];
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif //LICOM3_KOKKOS_SRC_FORTRAN_FORC_MOD_H_
