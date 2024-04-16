#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_PCONST_MOD_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_PCONST_MOD_H_

#include "def-undef.h"
#include "cpp_param_mod.h"

namespace CppPconstMod {
using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KMP1;
using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;

using CppParamMod::JMT_GLOBAL;
using CppParamMod::IMT_GLOBAL;

using CppParamMod::S_IMT;
using CppParamMod::S_JMT;
using CppParamMod::NTRA;
extern int (&j_global)[JMT];
extern int (&i_global)[IMT]; 

extern double *odz_pt;
extern int &ix;
extern int &iy;

extern double (&vit)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&viv)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double (&ahv_back)[JMT_GLOBAL];

extern int (&na)[MAX_BLOCKS_CLINIC][JMT][IMT];

extern double &dfricmx;
extern double &dwndmix;
extern double &pr_number;

#ifdef BCKMEX
extern double &diff_back_eq;
extern double &diff_back_coef;
extern double &diff_back_coef_max;
#endif // BCKMEX

#if (defined NETCDF) || (defined ALL)
extern float (&lon)[IMT_GLOBAL];
extern float (&lat)[JMT_GLOBAL];
extern float (&lon_o)[JMT_GLOBAL][IMT_GLOBAL];
extern float (&lat_o)[JMT_GLOBAL][IMT_GLOBAL];
extern float (&ulon_o)[JMT_GLOBAL][IMT_GLOBAL];
extern float (&ulat_o)[JMT_GLOBAL][IMT_GLOBAL];
extern float (&lev)[KM];
extern float (&lev1)[KM + 1];
#endif // NETCDF || ALL

extern double *s_lon;
extern double *s_lat;

extern double (&zkt)[KM];
extern double (&dzp)[KM];
extern double (&odzp)[KM];
extern double (&odzt)[KM];

extern double (&zkp)[KMP1];

extern double (&ebea)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&ebeb)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&ebla)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&eblb)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&epea)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&epeb)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&epla)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&eplb)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&rrd1)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&rrd2)[MAX_BLOCKS_CLINIC][JMT][IMT];

#if ( defined SMAG)

#endif // SMAG

extern double (&ohbt)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&ohbu)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&dzph)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&hbx)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&hby)[MAX_BLOCKS_CLINIC][JMT][IMT];
extern double (&snlat)[MAX_BLOCKS_CLINIC][JMT][IMT];

// extern int *&i_start;
// extern int *&j_start;

extern double (&to)[KM];
extern double (&so)[KM];
extern double (&c)[9][KM];
extern double (&po)[KM];

extern int &isop;

extern double &amv;
extern double &ahv;
extern double &ahice;

extern double (&akmu)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&akmt)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&akt)[MAX_BLOCKS_CLINIC][NTRA][KM][JMT][IMT];

extern double (&am)[JMT];
extern double (&ah)[JMT];

extern double (&am3)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&ah3)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double (&amx)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&amy)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double &gamma;

#if ( defined SMAG)
extern double &d0;
extern double &cp;
extern double &c0f;
extern double &tbice;
extern double &od0;
extern double &sag;
extern double &cag;
extern double &od0cp;
extern double &asea;
extern double &vsea;
extern double &afb1;
extern double &afb2;
extern double &afc1;
extern double &afc2;
extern double &aft1;
extern double &aft2;
extern double &karman;
extern double &rr;
#else
extern double &d0;
extern double &cp;
extern double &c0f;
extern double &tbice;
extern double &od0;
extern double &sag;
extern double &cag;
extern double &od0cp;
extern double &asea;
extern double &vsea;
extern double &afb1;
extern double &afb2;
extern double &afc1;
extern double &afc2;
extern double &aft1;
extern double &aft2;
extern double &karman;
extern double &rr;
#endif // SMAG

extern bool &diag_msf;
extern bool &diag_bsf;
extern bool &diag_mth;
extern bool &diag_budget;
extern bool &test_input;

extern char (&abmon)[12][3];
extern char (&abmon1)[12][3];

extern char (&out_dir)[80];

extern double &dtb;
extern double &dtc;
extern double &dts;
extern double &dtb2;
extern double &dtc2;
extern double &onbb;
extern double &onbc;
extern double &oncc;

extern int &nbb;
extern int &ncc;
extern int &nss;
extern int &isb;
extern int &isc;
extern int &ist;
extern int &month;

extern int (&nmonth)[12];
extern int (&nnmonth)[12];
extern int &number_day;
extern int &number_month;

extern int &number;
extern int &nstart; 
extern int &iy0;
extern int &iyfm;
extern int &mon0;
extern int &mend;
extern int &imd;
extern int &iday;
extern int &ii;
extern int &jj;
extern int &io_hist;
extern int &io_rest;
extern int &rest_freq;
extern int &hist_freq;
extern int &refdate;
extern int &boundary_restore;

extern int &klv;
extern int &kvt;

extern char (&adv_momentum)[80];
extern char (&adv_tracer)[80];

extern int &curr_ymd_licom;
extern int &curr_ymd_cpl;

extern int &yy_licom;
extern int &mm_licom;
extern int &dd_licom;
extern int &tod_licom;

extern int &yy_cpl;
extern int &mm_cpl;
extern int &dd_cpl;
extern int &tod_cpl;

extern int &imonth;
extern int &ihour;
extern int &iminute;
extern int &isecond;
extern int &tod;
extern int &ymd;
extern int &ymd_sync;
extern int &tod_sync;

extern int &ocn_cpl_dt;
extern int &licom_cpl_dt;

extern int &yearadd;
extern int &dayadd;

extern int &iyfmfnew;
extern int &first_step;
extern int &num_step_per_day;
extern int &num_output;
extern int &num_outputacc;
extern int &num_dailyoutput;
extern bool &dts_accum;
extern bool &daily_accum;
extern bool &dailybudget_accum;
extern bool &simple_assm;

#ifdef TIDEMIX
constexpr double LOCAL_MIXING_FRACTION = 0.33;
constexpr double MIXING_EF             = 0.20;
constexpr double SHELF_CUTOFF          = 1500.0;
constexpr double DECAY_SCALE           = 1000.0;
constexpr double BACK_TIDALMIXING      = 1.0e-6;
constexpr double MAX_TIDALMIXING       = 1.0e-2;
constexpr double MAX_WAVEDIS           = 0.1;
extern double (&fz_tide)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&ak_tide)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
  
extern double (&fztidal)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&richardson)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp3_tidal)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&ak_tide1)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

#endif // TIDEMIX

#ifdef CANUTOMIXOUT

extern double (&wp1_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp2_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp3_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp4_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp5_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp6_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp7_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp8_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp10_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp11_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp12_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wp13_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];

extern double (&wk1_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wk2_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wk3_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&wk4_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&fcor_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&fcort_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&alpha_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
extern double (&beta_canuto)[MAX_BLOCKS_CLINIC][KM][JMT][IMT];
#endif // CANUTOMIXOUT
} // namespace CppPconstMod

#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_PCONST_MOD_H_
