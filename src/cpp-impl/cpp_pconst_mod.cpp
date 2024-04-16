#include "../head/def-undef.h"

#include "../head/cpp_param_mod.h"
#include "../head/fortran_pconst_mod.h"

namespace CppPconstMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::S_IMT;
using CppParamMod::S_JMT;
using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;

double *s_lon = nullptr;
double *s_lat = nullptr;

double *odz_pt = nullptr;
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
int (&j_global)[JMT]                                    = pconst_mod_mp_j_global_; 
int (&i_global)[IMT]                                    = pconst_mod_mp_i_global_; 

int &ix                                                 = pconst_mod_mp_ix_;
int &iy                                                 = pconst_mod_mp_iy_;

double (&vit)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = pconst_mod_mp_vit_;
double (&viv)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = pconst_mod_mp_viv_;

// double (&ahv_back)[JMT_GLOBAL]                          = pconst_mod_mp_ahv_back_;

// int (&na)[MAX_BLOCKS_CLINIC][JMT][IMT]                  = pconst_mod_mp_na_;

double &dfricmx                                         = pconst_mod_mp_dfricmx_;
double &dwndmix                                         = pconst_mod_mp_dwndmix_;
double &pr_number                                       = pconst_mod_mp_pr_number_;

#ifdef BCKMEX
double &diff_back_eq                                    = pconst_mod_mp_diff_back_eq_;
double &diff_back_coef                                  = pconst_mod_mp_diff_back_coef_;
double &diff_back_coef_max                              = pconst_mod_mp_diff_back_coef_max_;
#endif // BCKMEX

#if (defined NETCDF) || (defined ALL)
float (&lon)[IMT_GLOBAL]                                = pconst_mod_mp_lon_;
float (&lat)[JMT_GLOBAL]                                = pconst_mod_mp_lat_;
float (&lon_o)[JMT_GLOBAL][IMT_GLOBAL]                  = pconst_mod_mp_lon_o_;
float (&lat_o)[JMT_GLOBAL][IMT_GLOBAL]                  = pconst_mod_mp_lat_o_;
float (&ulon_o)[JMT_GLOBAL][IMT_GLOBAL]                 = pconst_mod_mp_ulon_o_;
float (&ulat_o)[JMT_GLOBAL][IMT_GLOBAL]                 = pconst_mod_mp_ulat_o_;
float (&lev)[KM]                                        = pconst_mod_mp_lev_;
float (&lev1)[KM + 1]                                   = pconst_mod_mp_lev1_;
#endif // NETCDF || ALL

double (&zkt)[KM]                                       = pconst_mod_mp_zkt_;
double (&dzp)[KM]                                       = pconst_mod_mp_dzp_;
double (&odzp)[KM]                                      = pconst_mod_mp_odzp_;
double (&odzt)[KM]                                      = pconst_mod_mp_odzt_;

double (&zkp)[KMP1]                                     = pconst_mod_mp_zkp_;

double (&ebea)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_ebea_;
double (&ebeb)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_ebeb_;
double (&ebla)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_ebla_;
double (&eblb)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_eblb_;
double (&epea)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_epea_;
double (&epeb)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_epeb_;
double (&epla)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_epla_;
double (&eplb)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_eplb_;
double (&rrd1)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_rrd1_;
double (&rrd2)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_rrd2_;

#if ( defined SMAG)

#endif // SMAG

double (&ohbt)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_ohbt_;
double (&ohbu)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_ohbu_;
double (&dzph)[MAX_BLOCKS_CLINIC][JMT][IMT]             = pconst_mod_mp_dzph_;
double (&hbx)[MAX_BLOCKS_CLINIC][JMT][IMT]              = pconst_mod_mp_hbx_;
double (&hby)[MAX_BLOCKS_CLINIC][JMT][IMT]              = pconst_mod_mp_hby_;

double (&snlat)[MAX_BLOCKS_CLINIC][JMT][IMT]            = pconst_mod_mp_snlat_;

// int *(&i_start)                                         = pconst_mod_mp_i_start_;
// int *(&j_start)                                         = pconst_mod_mp_j_start_;

double (&to)[KM]                                        = pconst_mod_mp_to_;
double (&so)[KM]                                        = pconst_mod_mp_so_;
double (&c)[9][KM]                                      = pconst_mod_mp_c_;
double (&po)[KM]                                        = pconst_mod_mp_po_;

int &isop                                               = pconst_mod_mp_isop_;

double &amv                                             = pconst_mod_mp_amv_;
double &ahv                                             = pconst_mod_mp_ahv_;
double &ahice                                           = pconst_mod_mp_ahice_;

double (&akmu)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]         = pconst_mod_mp_akmu_;
double (&akmt)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]         = pconst_mod_mp_akmt_;
double (&akt)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = pconst_mod_mp_akt_;

double (&am)[JMT]                                       = pconst_mod_mp_am_;
double (&ah)[JMT]                                       = pconst_mod_mp_ah_;

double (&am3)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = pconst_mod_mp_am3_;
// double (&ah3)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = pconst_mod_mp_ah3_;

// double (&amx)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = pconst_mod_mp_amx_;
// double (&amy)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = pconst_mod_mp_amy_;

double &gamma                                           = pconst_mod_mp_gamma_;

#ifdef SMAG
double &d0                                              = pconst_mod_mp_d0_;
double &cp                                              = pconst_mod_mp_cp_;
double &c0f                                             = pconst_mod_mp_c0f_;
double &tbice                                           = pconst_mod_mp_tbice_;
double &od0                                             = pconst_mod_mp_od0_;
double &sag                                             = pconst_mod_mp_sag_;
double &cag                                             = pconst_mod_mp_cag_;
double &od0cp                                           = pconst_mod_mp_od0cp_;
double &asea                                            = pconst_mod_mp_asea_;
double &vsea                                            = pconst_mod_mp_vsea_;
double &afb1                                            = pconst_mod_mp_afb1_;
double &afb2                                            = pconst_mod_mp_afb2_;
double &afc1                                            = pconst_mod_mp_afc1_;
double &afc2                                            = pconst_mod_mp_afc2_;
double &aft1                                            = pconst_mod_mp_aft1_;
double &aft2                                            = pconst_mod_mp_aft2_;
double &karman                                          = pconst_mod_mp_karman_;
double &rr                                              = pconst_mod_mp_rr_;
#else // SMAG
double &d0                                              = pconst_mod_mp_d0_;
double &cp                                              = pconst_mod_mp_cp_;
double &c0f                                             = pconst_mod_mp_c0f_;
double &tbice                                           = pconst_mod_mp_tbice_;
double &od0                                             = pconst_mod_mp_od0_;
double &sag                                             = pconst_mod_mp_sag_;
double &cag                                             = pconst_mod_mp_cag_;
double &od0cp                                           = pconst_mod_mp_od0cp_;
double &asea                                            = pconst_mod_mp_asea_;
double &vsea                                            = pconst_mod_mp_vsea_;
double &afb1                                            = pconst_mod_mp_afb1_;
double &afb2                                            = pconst_mod_mp_afb2_;
double &afc1                                            = pconst_mod_mp_afc1_;
double &afc2                                            = pconst_mod_mp_afc2_;
double &aft1                                            = pconst_mod_mp_aft1_;
double &aft2                                            = pconst_mod_mp_aft2_;
#endif // SMAG

bool &diag_msf                                          = pconst_mod_mp_diag_msf_;
bool &diag_bsf                                          = pconst_mod_mp_diag_bsf_;
bool &diag_mth                                          = pconst_mod_mp_diag_mth_;
bool &diag_budget                                       = pconst_mod_mp_diag_budget_;
bool &test_input                                        = pconst_mod_mp_test_input_;

char (&abmon)[12][3]                                    = pconst_mod_mp_abmon_;
char (&abmon1)[12][3]                                   = pconst_mod_mp_abmon1_;

char (&out_dir)[80]                                     = pconst_mod_mp_out_dir_;

double &dtb                                             = pconst_mod_mp_dtb_;
double &dtc                                             = pconst_mod_mp_dtc_;
double &dts                                             = pconst_mod_mp_dts_;
double &dtb2                                            = pconst_mod_mp_dtb2_;
double &dtc2                                            = pconst_mod_mp_dtc2_;
double &onbb                                            = pconst_mod_mp_onbb_;
double &onbc                                            = pconst_mod_mp_onbc_;
double &oncc                                            = pconst_mod_mp_oncc_;

int &nbb                                                = pconst_mod_mp_nbb_;
int &ncc                                                = pconst_mod_mp_ncc_;
int &nss                                                = pconst_mod_mp_nss_;
int &isb                                                = pconst_mod_mp_isb_;
int &isc                                                = pconst_mod_mp_isc_;
int &ist                                                = pconst_mod_mp_ist_;
int &month                                              = pconst_mod_mp_month_;

int (&nmonth)[12]                                       = pconst_mod_mp_nmonth_;
int (&nnmonth)[12]                                      = pconst_mod_mp_nnmonth_;
int &number_day                                         = pconst_mod_mp_number_day_;
int &number_month                                       = pconst_mod_mp_number_month_;

int &number                                             = pconst_mod_mp_number_;
int &nstart                                             = pconst_mod_mp_nstart_;
int &iy0                                                = pconst_mod_mp_iy0_;
int &iyfm                                               = pconst_mod_mp_iyfm_;
int &mon0                                               = pconst_mod_mp_mon0_;
int &mend                                               = pconst_mod_mp_mend_;
int &imd                                                = pconst_mod_mp_imd_;
int &iday                                               = pconst_mod_mp_iday_;
int &ii                                                 = pconst_mod_mp_ii_;
int &jj                                                 = pconst_mod_mp_jj_;
int &io_hist                                            = pconst_mod_mp_io_hist_;
int &io_rest                                            = pconst_mod_mp_io_rest_;
int &rest_freq                                          = pconst_mod_mp_rest_freq_;
int &hist_freq                                          = pconst_mod_mp_hist_freq_;
int &refdate                                            = pconst_mod_mp_refdate_;
int &boundary_restore                                   = pconst_mod_mp_boundary_restore_;

int &klv                                                = pconst_mod_mp_klv_;
int &kvt                                                = pconst_mod_mp_kvt_;

char (&adv_momentum)[80]                                = pconst_mod_mp_adv_momentum_;
char (&adv_tracer)[80]                                  = pconst_mod_mp_adv_tracer_;

int &curr_ymd_licom                                     = pconst_mod_mp_curr_ymd_licom_;
int &curr_ymd_cpl                                       = pconst_mod_mp_curr_ymd_cpl_;

int &yy_licom                                           = pconst_mod_mp_yy_licom_;
int &mm_licom                                           = pconst_mod_mp_mm_licom_;
int &dd_licom                                           = pconst_mod_mp_dd_licom_;
int &tod_licom                                          = pconst_mod_mp_tod_licom_;

int &yy_cpl                                             = pconst_mod_mp_yy_cpl_;
int &mm_cpl                                             = pconst_mod_mp_mm_cpl_;
int &dd_cpl                                             = pconst_mod_mp_dd_cpl_;
int &tod_cpl                                            = pconst_mod_mp_tod_cpl_;

int &imonth                                             = pconst_mod_mp_imonth_;
int &ihour                                              = pconst_mod_mp_ihour_;
int &iminute                                            = pconst_mod_mp_iminute_;
int &isecond                                            = pconst_mod_mp_isecond_;
int &tod                                                = pconst_mod_mp_tod_;
int &ymd                                                = pconst_mod_mp_ymd_;
int &ymd_sync                                           = pconst_mod_mp_ymd_sync_;
int &tod_sync                                           = pconst_mod_mp_tod_sync_;

int &ocn_cpl_dt                                         = pconst_mod_mp_ocn_cpl_dt_;
int &licom_cpl_dt                                       = pconst_mod_mp_licom_cpl_dt_;

int &yearadd                                            = pconst_mod_mp_yearadd_;
int &dayadd                                             = pconst_mod_mp_dayadd_;

int &iyfmfnew                                           = pconst_mod_mp_iyfmfnew_;
int &first_step                                         = pconst_mod_mp_first_step_;
int &num_step_per_day                                   = pconst_mod_mp_num_step_per_day_;
int &num_output                                         = pconst_mod_mp_num_output_;
int &num_outputacc                                      = pconst_mod_mp_num_outputacc_;
int &num_dailyoutput                                    = pconst_mod_mp_num_dailyoutput_;
bool &dts_accum                                         = pconst_mod_mp_dts_accum_;
bool &daily_accum                                       = pconst_mod_mp_daily_accum_;
bool &dailybudget_accum                                 = pconst_mod_mp_dailybudget_accum_;
bool &simple_assm                                       = pconst_mod_mp_simple_assm_;

#ifdef TIDEMIX
double (&fz_tide)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]      = pconst_mod_mp_fz_tide_;
double (&ak_tide)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]      = pconst_mod_mp_ak_tide_;
  
double (&fztidal)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]      = pconst_mod_mp_fztidal_;
double (&richardson)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_richardson_;
double (&wp3_tidal)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]    = pconst_mod_mp_wp3_tidal_;
double (&ak_tide1)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]     = pconst_mod_mp_ak_tide1_;

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

double (&wp1_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wp1_canuto_;
double (&wp2_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wp2_canuto_;
double (&wp3_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wp3_canuto_;
double (&wp4_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wp4_canuto_;
double (&wp5_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wp5_canuto_;
double (&wp6_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wp6_canuto_;
double (&wp7_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wp7_canuto_;
double (&wp8_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wp8_canuto_;
double (&wp10_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = pconst_mod_mp_wp10_canuto_;
double (&wp11_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = pconst_mod_mp_wp11_canuto_;
double (&wp12_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = pconst_mod_mp_wp12_canuto_;
double (&wp13_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = pconst_mod_mp_wp13_canuto_;

double (&wk1_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wk1_canuto_;
double (&wk2_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wk2_canuto_;
double (&wk3_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wk3_canuto_;
double (&wk4_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = pconst_mod_mp_wk4_canuto_;
double (&fcor_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = pconst_mod_mp_fcor_canuto_;
double (&fcort_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = pconst_mod_mp_fcort_canuto_;
double (&alpha_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = pconst_mod_mp_alpha_canuto_;
double (&beta_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = pconst_mod_mp_beta_canuto_;

#endif // CANUTOMIXOUT
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
int (&j_global)[JMT]                                    = __pconst_mod_MOD_j_global; 
int (&i_global)[IMT]                                    = __pconst_mod_MOD_i_global; 

int &ix                                                 = __pconst_mod_MOD_ix;
int &iy                                                 = __pconst_mod_MOD_iy;

double (&vit)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = __pconst_mod_MOD_vit;
double (&viv)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = __pconst_mod_MOD_viv;

// double (&ahv_back)[JMT_GLOBAL]                          = __pconst_mod_MOD_ahv_back;

// int (&na)[MAX_BLOCKS_CLINIC][JMT][IMT]                  = __pconst_mod_MOD_na;

double &dfricmx                                         = __pconst_mod_MOD_dfricmx;
double &dwndmix                                         = __pconst_mod_MOD_dwndmix;
double &pr_number                                       = __pconst_mod_MOD_pr_number;

#ifdef BCKMEX
double &diff_back_eq                                    = __pconst_mod_MOD_diff_back_eq;
double &diff_back_coef                                  = __pconst_mod_MOD_diff_back_coef;
double &diff_back_coef_max                              = __pconst_mod_MOD_diff_back_coef_max;
#endif // BCKMEX

#if (defined NETCDF) || (defined ALL)
float (&lon)[IMT_GLOBAL]                                = __pconst_mod_MOD_lon;
float (&lat)[JMT_GLOBAL]                                = __pconst_mod_MOD_lat;
float (&lon_o)[JMT_GLOBAL][IMT_GLOBAL]                  = __pconst_mod_MOD_lon_o;
float (&lat_o)[JMT_GLOBAL][IMT_GLOBAL]                  = __pconst_mod_MOD_lat_o;
float (&ulon_o)[JMT_GLOBAL][IMT_GLOBAL]                 = __pconst_mod_MOD_ulon_o;
float (&ulat_o)[JMT_GLOBAL][IMT_GLOBAL]                 = __pconst_mod_MOD_ulat_o;
float (&lev)[KM]                                        = __pconst_mod_MOD_lev;
float (&lev1)[KM + 1]                                   = __pconst_mod_MOD_lev1;
#endif // NETCDF || ALL

// double (&s_lon)[S_JMT][S_IMT]                           = __pconst_mod_MOD_s_lon;
// double (&s_lat)[S_JMT][S_IMT]                           = __pconst_mod_MOD_s_lat;

double (&zkt)[KM]                                       = __pconst_mod_MOD_zkt;
double (&dzp)[KM]                                       = __pconst_mod_MOD_dzp;
double (&odzp)[KM]                                      = __pconst_mod_MOD_odzp;
double (&odzt)[KM]                                      = __pconst_mod_MOD_odzt;

double (&zkp)[KMP1]                                     = __pconst_mod_MOD_zkp;

double (&ebea)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_ebea;
double (&ebeb)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_ebeb;
double (&ebla)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_ebla;
double (&eblb)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_eblb;
double (&epea)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_epea;
double (&epeb)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_epeb;
double (&epla)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_epla;
double (&eplb)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_eplb;
double (&rrd1)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_rrd1;
double (&rrd2)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_rrd2;

#if ( defined SMAG)

#endif // SMAG

double (&ohbt)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_ohbt;
double (&ohbu)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_ohbu;
double (&dzph)[MAX_BLOCKS_CLINIC][JMT][IMT]             = __pconst_mod_MOD_dzph;
double (&hbx)[MAX_BLOCKS_CLINIC][JMT][IMT]              = __pconst_mod_MOD_hbx;
double (&hby)[MAX_BLOCKS_CLINIC][JMT][IMT]              = __pconst_mod_MOD_hby;

double (&snlat)[MAX_BLOCKS_CLINIC][JMT][IMT]            = __pconst_mod_MOD_snlat;

// int *(&i_start)                                         = __pconst_mod_MOD_i_start;
// int *(&j_start)                                         = __pconst_mod_MOD_j_start;

double (&to)[KM]                                        = __pconst_mod_MOD_to;
double (&so)[KM]                                        = __pconst_mod_MOD_so;
double (&c)[9][KM]                                      = __pconst_mod_MOD_c;
double (&po)[KM]                                        = __pconst_mod_MOD_po;

int &isop                                               = __pconst_mod_MOD_isop;

double &amv                                             = __pconst_mod_MOD_amv;
double &ahv                                             = __pconst_mod_MOD_ahv;
double &ahice                                           = __pconst_mod_MOD_ahice;

double (&akmu)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]         = __pconst_mod_MOD_akmu;
double (&akmt)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]         = __pconst_mod_MOD_akmt;
double (&akt)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT]    = __pconst_mod_MOD_akt;

double (&am)[JMT]                                       = __pconst_mod_MOD_am;
double (&ah)[JMT]                                       = __pconst_mod_MOD_ah;

double (&am3)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = __pconst_mod_MOD_am3;
// double (&ah3)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = __pconst_mod_MOD_ah3;

// double (&amx)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = __pconst_mod_MOD_amx;
// double (&amy)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]          = __pconst_mod_MOD_amy;

double &gamma                                           = __pconst_mod_MOD_gamma;

#ifdef SMAG
double &d0                                              = __pconst_mod_MOD_d0;
double &cp                                              = __pconst_mod_MOD_cp;
double &c0f                                             = __pconst_mod_MOD_c0f;
double &tbice                                           = __pconst_mod_MOD_tbice;
double &od0                                             = __pconst_mod_MOD_od0;
double &sag                                             = __pconst_mod_MOD_sag;
double &cag                                             = __pconst_mod_MOD_cag;
double &od0cp                                           = __pconst_mod_MOD_od0cp;
double &asea                                            = __pconst_mod_MOD_asea;
double &vsea                                            = __pconst_mod_MOD_vsea;
double &afb1                                            = __pconst_mod_MOD_afb1;
double &afb2                                            = __pconst_mod_MOD_afb2;
double &afc1                                            = __pconst_mod_MOD_afc1;
double &afc2                                            = __pconst_mod_MOD_afc2;
double &aft1                                            = __pconst_mod_MOD_aft1;
double &aft2                                            = __pconst_mod_MOD_aft2;
double &karman                                          = __pconst_mod_MOD_karman;
double &rr                                              = __pconst_mod_MOD_rr;
#else // SMAG
double &d0                                              = __pconst_mod_MOD_d0;
double &cp                                              = __pconst_mod_MOD_cp;
double &c0f                                             = __pconst_mod_MOD_c0f;
double &tbice                                           = __pconst_mod_MOD_tbice;
double &od0                                             = __pconst_mod_MOD_od0;
double &sag                                             = __pconst_mod_MOD_sag;
double &cag                                             = __pconst_mod_MOD_cag;
double &od0cp                                           = __pconst_mod_MOD_od0cp;
double &asea                                            = __pconst_mod_MOD_asea;
double &vsea                                            = __pconst_mod_MOD_vsea;
double &afb1                                            = __pconst_mod_MOD_afb1;
double &afb2                                            = __pconst_mod_MOD_afb2;
double &afc1                                            = __pconst_mod_MOD_afc1;
double &afc2                                            = __pconst_mod_MOD_afc2;
double &aft1                                            = __pconst_mod_MOD_aft1;
double &aft2                                            = __pconst_mod_MOD_aft2;
#endif // SMAG

bool &diag_msf                                          = __pconst_mod_MOD_diag_msf;
bool &diag_bsf                                          = __pconst_mod_MOD_diag_bsf;
bool &diag_mth                                          = __pconst_mod_MOD_diag_mth;
bool &diag_budget                                       = __pconst_mod_MOD_diag_budget;
bool &test_input                                        = __pconst_mod_MOD_test_input;

char (&abmon)[12][3]                                    = __pconst_mod_MOD_abmon;
char (&abmon1)[12][3]                                   = __pconst_mod_MOD_abmon1;

char (&out_dir)[80]                                     = __pconst_mod_MOD_out_dir;

double &dtb                                             = __pconst_mod_MOD_dtb;
double &dtc                                             = __pconst_mod_MOD_dtc;
double &dts                                             = __pconst_mod_MOD_dts;
double &dtb2                                            = __pconst_mod_MOD_dtb2;
double &dtc2                                            = __pconst_mod_MOD_dtc2;
double &onbb                                            = __pconst_mod_MOD_onbb;
double &onbc                                            = __pconst_mod_MOD_onbc;
double &oncc                                            = __pconst_mod_MOD_oncc;

int &nbb                                                = __pconst_mod_MOD_nbb;
int &ncc                                                = __pconst_mod_MOD_ncc;
int &nss                                                = __pconst_mod_MOD_nss;
int &isb                                                = __pconst_mod_MOD_isb;
int &isc                                                = __pconst_mod_MOD_isc;
int &ist                                                = __pconst_mod_MOD_ist;
int &month                                              = __pconst_mod_MOD_month;

int (&nmonth)[12]                                       = __pconst_mod_MOD_nmonth;
int (&nnmonth)[12]                                      = __pconst_mod_MOD_nnmonth;
int &number_day                                         = __pconst_mod_MOD_number_day;
int &number_month                                       = __pconst_mod_MOD_number_month;

int &number                                             = __pconst_mod_MOD_number;
int &nstart                                             = __pconst_mod_MOD_nstart;
int &iy0                                                = __pconst_mod_MOD_iy0;
int &iyfm                                               = __pconst_mod_MOD_iyfm;
int &mon0                                               = __pconst_mod_MOD_mon0;
int &mend                                               = __pconst_mod_MOD_mend;
int &imd                                                = __pconst_mod_MOD_imd;
int &iday                                               = __pconst_mod_MOD_iday;
int &ii                                                 = __pconst_mod_MOD_ii;
int &jj                                                 = __pconst_mod_MOD_jj;
int &io_hist                                            = __pconst_mod_MOD_io_hist;
int &io_rest                                            = __pconst_mod_MOD_io_rest;
int &rest_freq                                          = __pconst_mod_MOD_rest_freq;
int &hist_freq                                          = __pconst_mod_MOD_hist_freq;
int &refdate                                            = __pconst_mod_MOD_refdate;
int &boundary_restore                                   = __pconst_mod_MOD_boundary_restore;

int &klv                                                = __pconst_mod_MOD_klv;
int &kvt                                                = __pconst_mod_MOD_kvt;

char (&adv_momentum)[80]                                = __pconst_mod_MOD_adv_momentum;
char (&adv_tracer)[80]                                  = __pconst_mod_MOD_adv_tracer;

int &curr_ymd_licom                                     = __pconst_mod_MOD_curr_ymd_licom;
int &curr_ymd_cpl                                       = __pconst_mod_MOD_curr_ymd_cpl;

int &yy_licom                                           = __pconst_mod_MOD_yy_licom;
int &mm_licom                                           = __pconst_mod_MOD_mm_licom;
int &dd_licom                                           = __pconst_mod_MOD_dd_licom;
int &tod_licom                                          = __pconst_mod_MOD_tod_licom;

int &yy_cpl                                             = __pconst_mod_MOD_yy_cpl;
int &mm_cpl                                             = __pconst_mod_MOD_mm_cpl;
int &dd_cpl                                             = __pconst_mod_MOD_dd_cpl;
int &tod_cpl                                            = __pconst_mod_MOD_tod_cpl;

int &imonth                                             = __pconst_mod_MOD_imonth;
int &ihour                                              = __pconst_mod_MOD_ihour;
int &iminute                                            = __pconst_mod_MOD_iminute;
int &isecond                                            = __pconst_mod_MOD_isecond;
int &tod                                                = __pconst_mod_MOD_tod;
int &ymd                                                = __pconst_mod_MOD_ymd;
int &ymd_sync                                           = __pconst_mod_MOD_ymd_sync;
int &tod_sync                                           = __pconst_mod_MOD_tod_sync;

int &ocn_cpl_dt                                         = __pconst_mod_MOD_ocn_cpl_dt;
int &licom_cpl_dt                                       = __pconst_mod_MOD_licom_cpl_dt;

int &yearadd                                            = __pconst_mod_MOD_yearadd;
int &dayadd                                             = __pconst_mod_MOD_dayadd;

int &iyfmfnew                                           = __pconst_mod_MOD_iyfmfnew;
int &first_step                                         = __pconst_mod_MOD_first_step;
int &num_step_per_day                                   = __pconst_mod_MOD_num_step_per_day;
int &num_output                                         = __pconst_mod_MOD_num_output;
int &num_outputacc                                      = __pconst_mod_MOD_num_outputacc;
int &num_dailyoutput                                    = __pconst_mod_MOD_num_dailyoutput;
bool &dts_accum                                         = __pconst_mod_MOD_dts_accum;
bool &daily_accum                                       = __pconst_mod_MOD_daily_accum;
bool &dailybudget_accum                                 = __pconst_mod_MOD_dailybudget_accum;
bool &simple_assm                                       = __pconst_mod_MOD_simple_assm;

#ifdef TIDEMIX
double (&fz_tide)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]      = __pconst_mod_MOD_fz_tide;
double (&ak_tide)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]      = __pconst_mod_MOD_ak_tide;
  
double (&fztidal)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]      = __pconst_mod_MOD_fztidal;
double (&richardson)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_richardson;
double (&wp3_tidal)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]    = __pconst_mod_MOD_wp3_tidal;
double (&ak_tide1)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]     = __pconst_mod_MOD_ak_tide1;

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

double (&wp1_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wp1_canuto;
double (&wp2_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wp2_canuto;
double (&wp3_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wp3_canuto;
double (&wp4_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wp4_canuto;
double (&wp5_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wp5_canuto;
double (&wp6_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wp6_canuto;
double (&wp7_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wp7_canuto;
double (&wp8_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wp8_canuto;
double (&wp10_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __pconst_mod_MOD_wp10_canuto;
double (&wp11_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __pconst_mod_MOD_wp11_canuto;
double (&wp12_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __pconst_mod_MOD_wp12_canuto;
double (&wp13_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __pconst_mod_MOD_wp13_canuto;

double (&wk1_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wk1_canuto;
double (&wk2_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wk2_canuto;
double (&wk3_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wk3_canuto;
double (&wk4_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]   = __pconst_mod_MOD_wk4_canuto;
double (&fcor_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __pconst_mod_MOD_fcor_canuto;
double (&fcort_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __pconst_mod_MOD_fcort_canuto;
double (&alpha_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT] = __pconst_mod_MOD_alpha_canuto;
double (&beta_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]  = __pconst_mod_MOD_beta_canuto;

#endif // CANUTOMIXOUT
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
} // namespace CppPconstMod
