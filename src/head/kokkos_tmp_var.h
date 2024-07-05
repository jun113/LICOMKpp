#include "def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS
#include "kokkos_config.hpp"
namespace KokkosTmpVar {

// READYT
extern ViewDouble3D *p_v_work1;
extern ViewDouble3D *p_v_work2;

extern ViewDouble4D *p_v_pp;
extern ViewDouble4D *p_v_ppa;
extern ViewDouble4D *p_v_ppb;
extern ViewDouble4D *p_v_ppc;

extern ViewDouble4D *p_v_alpha;
extern ViewDouble4D *p_v_beta;

// READYC
extern ViewDouble4D *p_v_wp12;
extern ViewDouble4D *p_v_wp13;

extern ViewDouble3D *p_v_wk1;
extern ViewDouble3D *p_v_wk2;
extern ViewDouble3D *p_v_wk3;
extern ViewDouble3D *p_v_wk4;
extern ViewDouble3D *p_v_wp1;
extern ViewDouble3D *p_v_wp2;
extern ViewDouble3D *p_v_wp3;
extern ViewDouble3D *p_v_wp4;
extern ViewDouble3D *p_v_wp5;
extern ViewDouble3D *p_v_wp6;
extern ViewDouble3D *p_v_wp7;
extern ViewDouble3D *p_v_wp8;
extern ViewDouble3D *p_v_zlev;
extern ViewDouble3D * p_v_ak_tide_mixing;
// extern ViewInt2D *p_v_mld_lev;
extern ViewDouble3D *p_v_Ri;
extern ViewDouble3D *p_v_Rrho;
extern ViewDouble3D *p_v_Gm;

#ifdef BCKMEX
ViewDouble3D *p_v_diff_back    = nullptr;
ViewDouble3D *p_v_diff_back_sh = nullptr;
ViewDouble3D *p_v_diff_back_nn = nullptr;
#endif // BCKMEX
extern ViewDouble4D *p_v_uv_ws_face;

#ifndef SMAG

extern ViewDouble2D *p_v_div_out;

#ifdef BIHAR
extern ViewDouble2D *p_v_curl;

extern ViewDouble2D *p_v_cc;
extern ViewDouble3D *p_v_d2uk;
extern ViewDouble3D *p_v_d2vk;
#endif // BIHAR
#endif // SMAG

// TRACER
extern ViewDouble4D *p_v_vtl_ori;

extern ViewDouble3D *p_v_adv_tt;

extern ViewDouble4D *p_v_at_00_max_min;

#ifdef BIHAR
extern ViewDouble3D *p_v_dt2k;
#endif // BIHAR

extern ViewInt1D *p_v_nn;
extern ViewDouble1D *p_v_xs;

extern ViewDouble4D *p_v_c_cnsew;
#ifdef ISO
#ifdef LDD97
extern ViewDouble4D *p_v_f1;
extern ViewDouble4D *p_v_f2;
#endif // LDD97
#endif // ISO

// BAROTR
extern ViewDouble2D *p_v_gradx;
extern ViewDouble2D *p_v_grady;

// BCLINC
extern ViewDouble4D *p_v_work_merged;

// POP Halo Update
extern ViewDouble1D *p_v_halo_buffer;
} // namespace KokkosTmpVar

#endif // LICOM_ENABLE_KOKKOS
