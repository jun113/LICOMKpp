#ifndef LICOM3_KOKKOS_SRC_HEAD_FORTRAN_PCONST_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_FORTRAN_PCONST_MOD_H_
#include "def-undef.h"
#include "cpp_param_mod.h"

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMP1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::S_IMT;
using CppParamMod::S_JMT;
using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;
using CppParamMod::NTRA;
// extern
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
extern int pconst_mod_mp_j_global_[JMT]; 
extern int pconst_mod_mp_i_global_[IMT]; 

extern int pconst_mod_mp_ix_;
extern int pconst_mod_mp_iy_;

extern double pconst_mod_mp_vit_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_viv_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double pconst_mod_mp_ahv_back_[JMT_GLOBAL];

extern int pconst_mod_mp_na_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double pconst_mod_mp_dfricmx_;
extern double pconst_mod_mp_dwndmix_;
extern double pconst_mod_mp_pr_number_;

#ifdef BCKMEX
extern double pconst_mod_mp_diff_back_eq_;
extern double pconst_mod_mp_diff_back_coef_;
extern double pconst_mod_mp_diff_back_coef_max_;
#endif // BCKMEX

#if (defined NETCDF) || (defined ALL)
extern float pconst_mod_mp_lon_[IMT_GLOBAL];
extern float pconst_mod_mp_lat_[JMT_GLOBAL];
extern float pconst_mod_mp_lon_o_[JMT_GLOBAL][IMT_GLOBAL];
extern float pconst_mod_mp_lat_o_[JMT_GLOBAL][IMT_GLOBAL];
extern float pconst_mod_mp_ulon_o_[JMT_GLOBAL][IMT_GLOBAL];
extern float pconst_mod_mp_ulat_o_[JMT_GLOBAL][IMT_GLOBAL];
extern float pconst_mod_mp_lev_[KM];
extern float pconst_mod_mp_lev1_[KM + 1];
#endif // NETCDF || ALL

// extern double pconst_mod_mp_s_lon_[S_JMT][S_IMT];
// extern double pconst_mod_mp_s_lat_[S_JMT][S_IMT];

extern double pconst_mod_mp_zkt_[KM];
extern double pconst_mod_mp_dzp_[KM];
extern double pconst_mod_mp_odzp_[KM];
extern double pconst_mod_mp_odzt_[KM];

extern double pconst_mod_mp_zkp_[KMP1];

extern double pconst_mod_mp_ebea_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_ebeb_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_ebla_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_eblb_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_epea_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_epeb_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_epla_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_eplb_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double pconst_mod_mp_rrd1_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_rrd2_[MAX_BLOCKS_CLINIC][JMT][IMT];

#if ( defined SMAG)

#endif // SMAG

extern double pconst_mod_mp_ohbt_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_ohbu_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_dzph_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_hbx_[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double pconst_mod_mp_hby_[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double pconst_mod_mp_snlat_[MAX_BLOCKS_CLINIC][JMT][IMT];

// extern int *pconst_mod_mp_i_start_;
// extern int *pconst_mod_mp_j_start_;

extern double pconst_mod_mp_to_[KM];
extern double pconst_mod_mp_so_[KM];
extern double pconst_mod_mp_c_[9][KM];
extern double pconst_mod_mp_po_[KM];

extern int pconst_mod_mp_isop_;

extern double pconst_mod_mp_amv_;
extern double pconst_mod_mp_ahv_;
extern double pconst_mod_mp_ahice_;

extern double pconst_mod_mp_akmu_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_akmt_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_akt_[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double pconst_mod_mp_am_[JMT];
extern double pconst_mod_mp_ah_[JMT];

extern double pconst_mod_mp_am3_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_ah3_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double pconst_mod_mp_amx_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_amy_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double pconst_mod_mp_gamma_;

#ifdef SMAG
extern double pconst_mod_mp_d0_;
extern double pconst_mod_mp_cp_;
extern double pconst_mod_mp_c0f_;
extern double pconst_mod_mp_tbice_;
extern double pconst_mod_mp_od0_;
extern double pconst_mod_mp_sag_;
extern double pconst_mod_mp_cag_;
extern double pconst_mod_mp_od0cp_;
extern double pconst_mod_mp_asea_;
extern double pconst_mod_mp_vsea_;
extern double pconst_mod_mp_afb1_;
extern double pconst_mod_mp_afb2_;
extern double pconst_mod_mp_afc1_;
extern double pconst_mod_mp_afc2_;
extern double pconst_mod_mp_aft1_;
extern double pconst_mod_mp_aft2_;
extern double pconst_mod_mp_karman_;
extern double pconst_mod_mp_rr_;
#else
extern double pconst_mod_mp_d0_;
extern double pconst_mod_mp_cp_;
extern double pconst_mod_mp_c0f_;
extern double pconst_mod_mp_tbice_;
extern double pconst_mod_mp_od0_;
extern double pconst_mod_mp_sag_;
extern double pconst_mod_mp_cag_;
extern double pconst_mod_mp_od0cp_;
extern double pconst_mod_mp_asea_;
extern double pconst_mod_mp_vsea_;
extern double pconst_mod_mp_afb1_;
extern double pconst_mod_mp_afb2_;
extern double pconst_mod_mp_afc1_;
extern double pconst_mod_mp_afc2_;
extern double pconst_mod_mp_aft1_;
extern double pconst_mod_mp_aft2_;
extern double pconst_mod_mp_karman_;
extern double pconst_mod_mp_rr_;
#endif // SMAG

extern bool pconst_mod_mp_diag_msf_;
extern bool pconst_mod_mp_diag_bsf_;
extern bool pconst_mod_mp_diag_mth_;
extern bool pconst_mod_mp_diag_budget_;
extern bool pconst_mod_mp_test_input_;

extern char pconst_mod_mp_abmon_[12][3];
extern char pconst_mod_mp_abmon1_[12][3];

extern char pconst_mod_mp_out_dir_[80];

extern double pconst_mod_mp_dtb_;
extern double pconst_mod_mp_dtc_;
extern double pconst_mod_mp_dts_;
extern double pconst_mod_mp_dtb2_;
extern double pconst_mod_mp_dtc2_;
extern double pconst_mod_mp_onbb_;
extern double pconst_mod_mp_onbc_;
extern double pconst_mod_mp_oncc_;

extern int pconst_mod_mp_nbb_;
extern int pconst_mod_mp_ncc_;
extern int pconst_mod_mp_nss_;
extern int pconst_mod_mp_isb_;
extern int pconst_mod_mp_isc_;
extern int pconst_mod_mp_ist_;
extern int pconst_mod_mp_month_;

extern int pconst_mod_mp_nmonth_[12];
extern int pconst_mod_mp_nnmonth_[12];
extern int pconst_mod_mp_number_day_;
extern int pconst_mod_mp_number_month_;

extern int pconst_mod_mp_number_;
extern int pconst_mod_mp_nstart_;
extern int pconst_mod_mp_iy0_;
extern int pconst_mod_mp_iyfm_;
extern int pconst_mod_mp_mon0_;
extern int pconst_mod_mp_mend_;
extern int pconst_mod_mp_imd_;
extern int pconst_mod_mp_iday_;
extern int pconst_mod_mp_ii_;
extern int pconst_mod_mp_jj_;
extern int pconst_mod_mp_io_hist_;
extern int pconst_mod_mp_io_rest_;
extern int pconst_mod_mp_rest_freq_;
extern int pconst_mod_mp_hist_freq_;
extern int pconst_mod_mp_refdate_;
extern int pconst_mod_mp_boundary_restore_;

extern int pconst_mod_mp_klv_;
extern int pconst_mod_mp_kvt_;

extern char pconst_mod_mp_adv_momentum_[80];
extern char pconst_mod_mp_adv_tracer_[80];

extern int pconst_mod_mp_curr_ymd_licom_;
extern int pconst_mod_mp_curr_ymd_cpl_;

extern int pconst_mod_mp_yy_licom_;
extern int pconst_mod_mp_mm_licom_;
extern int pconst_mod_mp_dd_licom_;
extern int pconst_mod_mp_tod_licom_;

extern int pconst_mod_mp_yy_cpl_;
extern int pconst_mod_mp_mm_cpl_;
extern int pconst_mod_mp_dd_cpl_;
extern int pconst_mod_mp_tod_cpl_;

extern int pconst_mod_mp_imonth_;
extern int pconst_mod_mp_ihour_;
extern int pconst_mod_mp_iminute_;
extern int pconst_mod_mp_isecond_;
extern int pconst_mod_mp_tod_;
extern int pconst_mod_mp_ymd_;
extern int pconst_mod_mp_ymd_sync_;
extern int pconst_mod_mp_tod_sync_;

extern int pconst_mod_mp_ocn_cpl_dt_;
extern int pconst_mod_mp_licom_cpl_dt_;

extern int pconst_mod_mp_yearadd_;
extern int pconst_mod_mp_dayadd_;

extern int pconst_mod_mp_iyfmfnew_;
extern int pconst_mod_mp_first_step_;
extern int pconst_mod_mp_num_step_per_day_;
extern int pconst_mod_mp_num_output_;
extern int pconst_mod_mp_num_outputacc_;
extern int pconst_mod_mp_num_dailyoutput_;
extern bool pconst_mod_mp_dts_accum_;
extern bool pconst_mod_mp_daily_accum_;
extern bool pconst_mod_mp_dailybudget_accum_;
extern bool pconst_mod_mp_simple_assm_;

#ifdef TIDEMIX
extern double pconst_mod_mp_fz_tide_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_ak_tide_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  
extern double pconst_mod_mp_fztidal_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_richardson_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp3_tidal_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_ak_tide1_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

extern double pconst_mod_mp_wp1_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp2_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp3_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp4_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp5_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp6_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp7_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp8_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp10_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp11_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp12_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wp13_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double pconst_mod_mp_wk1_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wk2_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wk3_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_wk4_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double pconst_mod_mp_fcor_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_fcort_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double pconst_mod_mp_alpha_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double pconst_mod_mp_beta_canuto_[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#endif // CANUTOMIXOUT
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
extern int __pconst_mod_MOD_j_global[JMT]; 
extern int __pconst_mod_MOD_i_global[IMT]; 

extern int __pconst_mod_MOD_ix;
extern int __pconst_mod_MOD_iy;

extern double __pconst_mod_MOD_vit[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_viv[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double __pconst_mod_MOD_ahv_back[JMT_GLOBAL];

extern int __pconst_mod_MOD_na[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __pconst_mod_MOD_dfricmx;
extern double __pconst_mod_MOD_dwndmix;
extern double __pconst_mod_MOD_pr_number;

#ifdef BCKMEX
extern double __pconst_mod_MOD_diff_back_eq;
extern double __pconst_mod_MOD_diff_back_coef;
extern double __pconst_mod_MOD_diff_back_coef_max;
#endif // BCKMEX

#if (defined NETCDF) || (defined ALL)
extern float __pconst_mod_MOD_lon[IMT_GLOBAL];
extern float __pconst_mod_MOD_lat[JMT_GLOBAL];
extern float __pconst_mod_MOD_lon_o[JMT_GLOBAL][IMT_GLOBAL];
extern float __pconst_mod_MOD_lat_o[JMT_GLOBAL][IMT_GLOBAL];
extern float __pconst_mod_MOD_ulon_o[JMT_GLOBAL][IMT_GLOBAL];
extern float __pconst_mod_MOD_ulat_o[JMT_GLOBAL][IMT_GLOBAL];
extern float __pconst_mod_MOD_lev[KM];
extern float __pconst_mod_MOD_lev1[KM + 1];
#endif // NETCDF || ALL

// extern double __pconst_mod_MOD_s_lon[S_JMT][S_IMT];
// extern double __pconst_mod_MOD_s_lat[S_JMT][S_IMT];

extern double __pconst_mod_MOD_zkt[KM];
extern double __pconst_mod_MOD_dzp[KM];
extern double __pconst_mod_MOD_odzp[KM];
extern double __pconst_mod_MOD_odzt[KM];

extern double __pconst_mod_MOD_zkp[KMP1];

extern double __pconst_mod_MOD_ebea[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_ebeb[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_ebla[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_eblb[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_epea[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_epeb[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_epla[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_eplb[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __pconst_mod_MOD_rrd1[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_rrd2[MAX_BLOCKS_CLINIC][JMT][IMT];

#if ( defined SMAG)

#endif // SMAG

extern double __pconst_mod_MOD_ohbt[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_ohbu[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_dzph[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_hbx[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double __pconst_mod_MOD_hby[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double __pconst_mod_MOD_snlat[MAX_BLOCKS_CLINIC][JMT][IMT];

// extern int *__pconst_mod_MOD_i_start;
// extern int *__pconst_mod_MOD_j_start;

extern double __pconst_mod_MOD_to[KM];
extern double __pconst_mod_MOD_so[KM];
extern double __pconst_mod_MOD_c[9][KM];
extern double __pconst_mod_MOD_po[KM];

extern int __pconst_mod_MOD_isop;

extern double __pconst_mod_MOD_amv;
extern double __pconst_mod_MOD_ahv;
extern double __pconst_mod_MOD_ahice;

extern double __pconst_mod_MOD_akmu[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_akmt[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_akt[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double __pconst_mod_MOD_am[JMT];
extern double __pconst_mod_MOD_ah[JMT];

extern double __pconst_mod_MOD_am3[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_ah3[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double __pconst_mod_MOD_amx[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_amy[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double __pconst_mod_MOD_gamma;

#ifdef SMAG
extern double __pconst_mod_MOD_d0;
extern double __pconst_mod_MOD_cp;
extern double __pconst_mod_MOD_c0f;
extern double __pconst_mod_MOD_tbice;
extern double __pconst_mod_MOD_od0;
extern double __pconst_mod_MOD_sag;
extern double __pconst_mod_MOD_cag;
extern double __pconst_mod_MOD_od0cp;
extern double __pconst_mod_MOD_asea;
extern double __pconst_mod_MOD_vsea;
extern double __pconst_mod_MOD_afb1;
extern double __pconst_mod_MOD_afb2;
extern double __pconst_mod_MOD_afc1;
extern double __pconst_mod_MOD_afc2;
extern double __pconst_mod_MOD_aft1;
extern double __pconst_mod_MOD_aft2;
extern double __pconst_mod_MOD_karman;
extern double __pconst_mod_MOD_rr;
#else
extern double __pconst_mod_MOD_d0;
extern double __pconst_mod_MOD_cp;
extern double __pconst_mod_MOD_c0f;
extern double __pconst_mod_MOD_tbice;
extern double __pconst_mod_MOD_od0;
extern double __pconst_mod_MOD_sag;
extern double __pconst_mod_MOD_cag;
extern double __pconst_mod_MOD_od0cp;
extern double __pconst_mod_MOD_asea;
extern double __pconst_mod_MOD_vsea;
extern double __pconst_mod_MOD_afb1;
extern double __pconst_mod_MOD_afb2;
extern double __pconst_mod_MOD_afc1;
extern double __pconst_mod_MOD_afc2;
extern double __pconst_mod_MOD_aft1;
extern double __pconst_mod_MOD_aft2;
extern double __pconst_mod_MOD_karman;
extern double __pconst_mod_MOD_rr;
#endif // SMAG

extern bool __pconst_mod_MOD_diag_msf;
extern bool __pconst_mod_MOD_diag_bsf;
extern bool __pconst_mod_MOD_diag_mth;
extern bool __pconst_mod_MOD_diag_budget;
extern bool __pconst_mod_MOD_test_input;

extern char __pconst_mod_MOD_abmon[12][3];
extern char __pconst_mod_MOD_abmon1[12][3];

extern char __pconst_mod_MOD_out_dir[80];

extern double __pconst_mod_MOD_dtb;
extern double __pconst_mod_MOD_dtc;
extern double __pconst_mod_MOD_dts;
extern double __pconst_mod_MOD_dtb2;
extern double __pconst_mod_MOD_dtc2;
extern double __pconst_mod_MOD_onbb;
extern double __pconst_mod_MOD_onbc;
extern double __pconst_mod_MOD_oncc;

extern int __pconst_mod_MOD_nbb;
extern int __pconst_mod_MOD_ncc;
extern int __pconst_mod_MOD_nss;
extern int __pconst_mod_MOD_isb;
extern int __pconst_mod_MOD_isc;
extern int __pconst_mod_MOD_ist;
extern int __pconst_mod_MOD_month;

extern int __pconst_mod_MOD_nmonth[12];
extern int __pconst_mod_MOD_nnmonth[12];
extern int __pconst_mod_MOD_number_day;
extern int __pconst_mod_MOD_number_month;

extern int __pconst_mod_MOD_number;
extern int __pconst_mod_MOD_nstart;
extern int __pconst_mod_MOD_iy0;
extern int __pconst_mod_MOD_iyfm;
extern int __pconst_mod_MOD_mon0;
extern int __pconst_mod_MOD_mend;
extern int __pconst_mod_MOD_imd;
extern int __pconst_mod_MOD_iday;
extern int __pconst_mod_MOD_ii;
extern int __pconst_mod_MOD_jj;
extern int __pconst_mod_MOD_io_hist;
extern int __pconst_mod_MOD_io_rest;
extern int __pconst_mod_MOD_rest_freq;
extern int __pconst_mod_MOD_hist_freq;
extern int __pconst_mod_MOD_refdate;
extern int __pconst_mod_MOD_boundary_restore;

extern int __pconst_mod_MOD_klv;
extern int __pconst_mod_MOD_kvt;

extern char __pconst_mod_MOD_adv_momentum[80];
extern char __pconst_mod_MOD_adv_tracer[80];

extern int __pconst_mod_MOD_curr_ymd_licom;
extern int __pconst_mod_MOD_curr_ymd_cpl;

extern int __pconst_mod_MOD_yy_licom;
extern int __pconst_mod_MOD_mm_licom;
extern int __pconst_mod_MOD_dd_licom;
extern int __pconst_mod_MOD_tod_licom;

extern int __pconst_mod_MOD_yy_cpl;
extern int __pconst_mod_MOD_mm_cpl;
extern int __pconst_mod_MOD_dd_cpl;
extern int __pconst_mod_MOD_tod_cpl;

extern int __pconst_mod_MOD_imonth;
extern int __pconst_mod_MOD_ihour;
extern int __pconst_mod_MOD_iminute;
extern int __pconst_mod_MOD_isecond;
extern int __pconst_mod_MOD_tod;
extern int __pconst_mod_MOD_ymd;
extern int __pconst_mod_MOD_ymd_sync;
extern int __pconst_mod_MOD_tod_sync;

extern int __pconst_mod_MOD_ocn_cpl_dt;
extern int __pconst_mod_MOD_licom_cpl_dt;

extern int __pconst_mod_MOD_yearadd;
extern int __pconst_mod_MOD_dayadd;

extern int __pconst_mod_MOD_iyfmfnew;
extern int __pconst_mod_MOD_first_step;
extern int __pconst_mod_MOD_num_step_per_day;
extern int __pconst_mod_MOD_num_output;
extern int __pconst_mod_MOD_num_outputacc;
extern int __pconst_mod_MOD_num_dailyoutput;
extern bool __pconst_mod_MOD_dts_accum;
extern bool __pconst_mod_MOD_daily_accum;
extern bool __pconst_mod_MOD_dailybudget_accum;
extern bool __pconst_mod_MOD_simple_assm;

#ifdef TIDEMIX
extern double __pconst_mod_MOD_fz_tide[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_ak_tide[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  
extern double __pconst_mod_MOD_fztidal[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_richardson[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp3_tidal[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_ak_tide1[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

extern double __pconst_mod_MOD_wp1_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp2_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp3_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp4_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp5_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp6_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp7_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp8_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp10_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp11_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp12_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wp13_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double __pconst_mod_MOD_wk1_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wk2_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wk3_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_wk4_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double __pconst_mod_MOD_fcor_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_fcort_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double __pconst_mod_MOD_alpha_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double __pconst_mod_MOD_beta_canuto[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#endif // CANUTOMIXOUT
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

#endif // LICOM3_KOKKOS_SRC_HEAD_FORTRAN_PCONST_MOD_H_
